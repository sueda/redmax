classdef (Abstract) Joint < handle
	%Joint A generic joint between two bodies
	%   A joint is defined between a parent and the child. The DOF, q, of
	%   the joint is the relative displacement of the child wrt the parent.
	
	%%
	properties
		name     % Optional name
		body     % Attached body
		parent   % Parent joint
		children % Children joints
		ndof     % Number of DOF
		q        % Position
		qdot     % Velocity
		qddot    % Acceleration
		qInit    % Initial position (for joint stiffness)
		tau      % Joint torque
		tauCon   % Constraint torque
		presc    % Prescribed motion constraint
		Kr       % Joint stiffness
		Dr       % Joint damping
		S        % Jacobian
		Sdot     % dS/dt
		I_j      % Inertia at the joint
		V        % Twist at parent joint
		Vdot     % Acceleration at parent joint
		E_pj     % Transform of this joint wrt parent joint
		E_pj0    % Transform when q is zero
		E_jp     % Transform of parent joint wrt this joint
		Ad_jp    % Adjoint of E_jp
		E_wj     % Transform of this joint wrt world
		contacts % Internal contacts
		next     % Forward recursive ordering
		prev     % Reverse recursive ordering
		idxR     % Reduced indices
		idxT     % Tangent matrix indices
	end
	
	%%
	properties (Access = protected)
		Q % Transformation matrix applied about the joint
	end
	
	%%
	properties (Access = private)
		% For Recursive Hybrid Dynamics
		Sqdot
		eta
		Ihat
		Bhat
		Psi
		Pi
		beta
		F
		% For J'*x product
		alpha
	end
	
	%%
	properties (Constant)
		qpOpt = optimoptions(@quadprog,'Display','off');
	end
	
	%%
	methods
		%%
		function this = Joint(parent,body,ndof)
			if isempty(parent)
				this.name = ['NULL-',body.name];
			else
				this.name = [parent.body.name,'-',body.name];
			end
			this.parent = parent;
			this.body = body;
			this.ndof = ndof;
			this.q = zeros(ndof,1);
			this.qdot = zeros(ndof,1);
			this.qddot = zeros(ndof,1);
			this.tau = zeros(ndof,1);
			this.tauCon = zeros(ndof,1);
			this.Kr = 0;
			this.Dr = 0;
			this.S = zeros(6,ndof);
			this.Sdot = zeros(6,ndof);
			this.I_j = eye(6);
			this.V = zeros(6,1);
			this.Vdot = zeros(6,1);
			this.contacts = {};
			body.joint = this;
			if ~isempty(parent)
				parent.children{end+1} = this;
			end
		end
		
		%%
		function setJointTransform(this,E)
			% Sets the transform of this joint wrt parent joint
			this.E_pj0 = E;
		end
		
		%%
		function setStiffness(this,stiffness)
			% Sets this joint's linear stiffness
			this.Kr = stiffness;
		end
		
		%%
		function setDamping(this,damping)
			% Sets this joint's linear velocity damping
			this.Dr = damping;
		end
		
		%%
		function order = getTraversalOrder(this,order)
			% Gets traversal order
			if nargin == 1
				order = [];
			end
			order(end+1) = length(order) + 1;
			for i = 1 : length(this.children)
				order = this.children{i}.getTraversalOrder(order);
			end
		end
		
		%%
		function countDofs(this)
			% Counts reduced DOFs
			nr = redmax.Scene.countR();
			this.idxR = nr + (1:this.ndof);
			nr = nr + this.ndof;
			redmax.Scene.countR(nr);
			this.body.countDofs();
			% Use the current q as the initial q
			this.qInit = this.q;
		end
		
		%%
		function computeInertia(this)
			m = this.body.I_i(4);
			R = this.body.E_ji(1:3,1:3);
			p = this.body.E_ji(1:3,4);
			pBrac = se3.brac(p);
			Ic = diag(this.body.I_i(1:3));
			this.I_j = [R*Ic*R'+m*(pBrac'*pBrac),m*pBrac;m*pBrac',m*eye(3)];
		end
		
		%%
		function reparam(this)
			% If necessary, reparameterizes q and qdot and then updates S
			% and Sdot appropriately. Accelerations, qddot, will not be
			% updated and will need to be computed by the integrator.
			this.reparam_();
			if ~isempty(this.next)
				this.next.reparam();
			end
		end
		
		%%
		function y = gatherDofs(this,y)
			% Gathers q and qdot into y
			nr = redmax.Scene.countR();
			if nargin == 1
				y = zeros(2*nr,1);
			end
			y(this.idxR) = this.q;
			y(nr + this.idxR) = this.qdot;
			if ~isempty(this.next)
				y = this.next.gatherDofs(y);
			end
		end
		
		%%
		function ydot = gatherDDofs(this,ydot)
			% Gathers qdot and qddot into ydot
			nr = redmax.Scene.countR();
			if nargin == 1
				ydot = zeros(2*nr,1);
			end
			ydot(this.idxR) = this.qdot;
			ydot(nr + this.idxR) = this.qddot;
			if ~isempty(this.next)
				ydot = this.next.gatherDDofs(ydot);
			end
		end
		
		%%
		function scatterDofs(this,y)
			% Scatters q and qdot from y
			nr = redmax.Scene.countR();
			this.q(1:this.ndof) = y(this.idxR);
			this.qdot(1:this.ndof) = y(nr + this.idxR);
			this.update(); % safe to call with forward traversal
			if ~isempty(this.next)
				this.next.scatterDofs(y);
			end
		end
		
		%%
		function scatterDDofs(this,ydot)
			% Scatters qdot and qddot from ydot
			nr = redmax.Scene.countR();
			this.qdot(1:this.ndof) = ydot(this.idxR);
			this.qddot(1:this.ndof) = ydot(nr + this.idxR);
			if ~isempty(this.next)
				this.next.scatterDDofs(ydot);
			end
		end
		
		%%
		function scatterTauCon(this,tauc)
			% Scatters constraint force
			if nargin == 1
				nr = redmax.Scene.countR();
				tauc = zeros(nr,1);
			end
			this.tauCon = tauc(this.idxR);
			if ~isempty(this.next)
				this.next.scatterTauCon(tauc);
			end
		end
		
		%%
		function update(this)
			% Updates this joint and the attached body
			this.update_();
			% Transforms and adjoints
			if isempty(this.E_pj0)
				this.E_pj = this.Q;
			else
				this.E_pj = this.E_pj0*this.Q;
			end
			this.E_jp = se3.inv(this.E_pj);
			this.Ad_jp = se3.Ad(this.E_jp);
			if isempty(this.parent)
				E_wp = eye(4);
			else
				E_wp = this.parent.E_wj;
			end
			this.E_wj = E_wp*this.E_pj;
			% Joint velocity
			if this.ndof == 0
				this.V = zeros(6,1);
			else
				this.V = this.S*this.qdot;
			end
			if ~isempty(this.parent) % Add parent velocity
				this.V = this.V + this.Ad_jp*this.parent.V;
			end
			this.body.update();
			if ~isempty(this.next)
				this.next.update();
			end
		end
		
		%%
		function rhdPass1(this)
			% Recursive hybrid dynamics pass 1
			
			% V
			this.Sqdot = this.S*this.qdot;
			this.V = this.Sqdot;
			if ~isempty(this.parent)
				this.V = this.V + this.Ad_jp*this.parent.V;
			end
			% eta
			this.eta = se3.ad(this.V)*this.Sqdot + this.Sdot*this.qdot;
			if ~isempty(this.next)
				this.next.rhdPass1();
			end
		end
		
		%%
		function rhdPass2(this,grav)
			% Recursive hybrid dynamics pass 2

			% Ihat
			this.Ihat = this.I_j;
			for k = 1 : length(this.children)
				child = this.children{k};
				this.Ihat = this.Ihat + child.Ad_jp'*child.Pi*child.Ad_jp;
			end
			% Bhat
			Fgrav = zeros(6,1);
			m = this.I_j(6,6);
			Fgrav(4:6) = this.body.E_wi(1:3,1:3)'*m*grav; % wrench in body space
			Fgrav = this.body.Ad_ij'*Fgrav;
			Fext = Fgrav;
			if ~isempty(this.body.wext_i)
				Fext = Fext + this.body.Ad_ij'*this.body.wext_i;
			end
			this.Bhat = -se3.ad(this.V)'*this.I_j*this.V - Fext;
			for k = 1 : length(this.children)
				child = this.children{k};
				this.Bhat = this.Bhat + child.Ad_jp'*child.beta;
			end
			if ~isempty(this.presc)
				% Pi
				this.qddot = this.presc.qddot;
				this.Pi = this.Ihat;
				this.beta = this.Bhat + this.Ihat*(this.eta + this.S*this.qddot);
			else
				% Psi
				this.Psi = inv(this.S'*this.Ihat*this.S);
				% Pi
				this.Pi = this.Ihat - this.Ihat*this.S*this.Psi*this.S'*this.Ihat;
				% beta
				tauTotal = this.tau - this.Kr*(this.q - this.qInit) - this.Dr*this.qdot + this.tauCon;
				this.beta = this.Bhat + this.Ihat*(this.eta + this.S*this.Psi*(tauTotal - this.S'*(this.Ihat*this.eta + this.Bhat)));
			end
			if ~isempty(this.prev)
				this.prev.rhdPass2(grav);
			end
		end
		
		%%
		function rhdPass3(this)
			% Recursive hybrid dynamics pass 3
			if ~isempty(this.presc)
				this.qddot = this.presc.qddot;
				% Vdot
				this.Vdot = zeros(6,1);
				if this.ndof > 0
					this.Vdot = this.S*this.qddot + this.eta;
				end
				if ~isempty(this.parent)
					this.Vdot = this.Vdot + this.Ad_jp*this.parent.Vdot;
				end
				% F
				this.F = this.Ihat*this.Vdot + this.Bhat;
				% tau
				this.tau = this.S'*this.F;
			else
				% qddot
				tmp = zeros(6,1);
				if ~isempty(this.parent)
					tmp = this.Ad_jp*this.parent.Vdot;
				end
				tauTotal = this.tau - this.Kr*(this.q - this.qInit) - this.Dr*this.qdot + this.tauCon;
				this.qddot = this.Psi*(tauTotal - this.S'*this.Ihat*(tmp + this.eta) - this.S'*this.Bhat);
				% Vdot
				this.Vdot = zeros(6,1);
				if this.ndof > 0
					this.Vdot = this.S*this.qddot + this.eta;
				end
				if ~isempty(this.parent)
					this.Vdot = this.Vdot + this.Ad_jp*this.parent.Vdot;
				end
				% F
				this.F = this.Ihat*this.Vdot + this.Bhat;
			end
			if ~isempty(this.next)
				this.next.rhdPass3();
			end
		end
		
		%%
		function computeMinvProdInit(this,h)
			% Initializes preconditioner solve.
			% Computes things that are the same for all RHS vectors
			% If the optional argument h is passed in, appends the damping
			% and stiffness values to the diagonal, for linearly implicit
			% Euler integration.
			%
			% WARNING: Ad_jp changes after taking a step, so inv(Mr) is
			% valid for the current step. Once scatter/update is called,
			% the inverse will not match.
			if nargin == 1
				h = 0;
			end
			% Ihat
			if h == 0
				% Just inertia
				this.Ihat = this.I_j;
			else
				% Add LHS terms for maximal damping and stiffness
				Ad_ij = this.body.Ad_ij;
				DKm = h*Ad_ij'*(this.body.Dmdiag + this.body.Ddiag - h*this.body.Kmdiag)*Ad_ij;
				this.Ihat = this.I_j + DKm;
			end
			for k = 1 : length(this.children)
				child = this.children{k};
				this.Ihat = this.Ihat + child.Ad_jp'*child.Pi*child.Ad_jp;
			end
			% Psi
			if h == 0
				% Just inertia
				this.Psi = inv(this.S'*this.Ihat*this.S);
			else
				% Add LHS terms for reduced damping and stiffness
				DKr = h*(this.Dr + h*this.Kr)*eye(this.ndof);
				this.Psi = inv(this.S'*this.Ihat*this.S + DKr);
			end
			% Pi
			this.Pi = this.Ihat - this.Ihat*this.S*this.Psi*this.S'*this.Ihat;
			if ~isempty(this.prev)
				this.prev.computeMinvProdInit(h);
			end
		end
		
		%%
		function computeMinvProd2(this,x)
			% Preconditioner solve pass 2
			x_i = x(this.idxR);
			% Bhat
			this.Bhat = zeros(6,1);
			for k = 1 : length(this.children)
				child = this.children{k};
				this.Bhat = this.Bhat + child.Ad_jp'*child.beta;
			end
			% beta
			this.beta = this.Bhat + this.Ihat*(this.S*this.Psi*(x_i - this.S'*this.Bhat));
			if ~isempty(this.prev)
				this.prev.computeMinvProd2(x);
			end
		end
		
		%%
		function y = computeMinvProd3(this,x,y)
			% Preconditioner solve pass 3
			if nargin == 2
				nr = redmax.Scene.countR();
				y = zeros(nr,1);
			end
			x_i = x(this.idxR);
			% qddot
			tmp = zeros(6,1);
			if ~isempty(this.parent)
				tmp = this.Ad_jp*this.parent.Vdot;
			end
			this.qddot = this.Psi*(x_i - this.S'*this.Ihat*tmp - this.S'*this.Bhat);
			y(this.idxR) = this.qddot;
			% Vdot
			this.Vdot = zeros(6,1);
			if this.ndof > 0
				this.Vdot = this.S*this.qddot;
			end
			if ~isempty(this.parent)
				this.Vdot = this.Vdot + this.Ad_jp*this.parent.Vdot;
			end
			if ~isempty(this.next)
				y = this.next.computeMinvProd3(x,y);
			end
		end
		
		%%
		function [fr,Kr] = computeForceStiffness(this,fr,Kr)
			% Computes joint stiffness force vector and matrix
			nr = redmax.Scene.countR();
			if nargin == 1
				fr = zeros(nr,1);
				Kr = zeros(nr);
			elseif nargin == 2
				Kr = zeros(nr);
			end
			if isempty(this.presc)
				rows = this.idxR;
				% We'll add the joint torque here rather than having a
				% separate function.
				fr(rows) = fr(rows) + this.tau - this.Kr*(this.q - this.qInit);
				Kr(rows,rows) = Kr(rows,rows) - this.Kr*eye(this.ndof);
			end
			if ~isempty(this.next)
				[fr,Kr] = this.next.computeForceStiffness(fr,Kr);
			end
		end
		
		%%
		function [fr,Dr] = computeForceDamping(this,fr,Dr)
			% Computes joint damping force vector and matrix
			nr = redmax.Scene.countR();
			if nargin == 1
				fr = zeros(nr,1);
				Dr = zeros(nr);
			elseif nargin == 2
				Dr = zeros(nr);
			end
			if isempty(this.presc)
				rows = this.idxR;
				fr(rows) = fr(rows) - this.Dr*this.qdot;
				Dr(rows,rows) = Dr(rows,rows) + this.Dr*eye(this.ndof);
			end
			if ~isempty(this.next)
				[fr,Dr] = this.next.computeForceDamping(fr,Dr);
			end
		end
		
		%%
		function [J,Jdot] = computeJacobian(this,J,Jdot)
			% Computes the redmax Jacobian
			nm = redmax.Scene.countM();
			nr = redmax.Scene.countR();
			if nargout == 1
				if nargin == 1
					J = zeros(nm,nr);
				end
				rows = this.body.idxM;
				cols = this.idxR;
				Ad_ij = this.body.Ad_ij;
				J(rows,cols) = Ad_ij*this.S;
				% Loop through all ancestors
				jointA = this.parent;
				while ~isempty(jointA)
					rowsP = this.parent.body.idxM;
					Ad_ip = this.body.Ad_ip;
					colsA = jointA.idxR;
					J(rows,colsA) = Ad_ip*J(rowsP,colsA);
					jointA = jointA.parent;
				end
				if ~isempty(this.next)
					J = this.next.computeJacobian(J);
				end
			else
				if nargin == 1
					J = zeros(nm,nr);
					Jdot = zeros(nm,nr);
				end
				rows = this.body.idxM;
				cols = this.idxR;
				Ad_ij = this.body.Ad_ij;
				J(rows,cols) = Ad_ij*this.S;
				Jdot(rows,cols) = Ad_ij*this.Sdot;
				% Loop through all ancestors
				jointA = this.parent;
				while ~isempty(jointA)
					rowsP = this.parent.body.idxM;
					Ad_ip = this.body.Ad_ip;
					Ad_iw = this.body.Ad_iw;
					Ad_wp = this.parent.body.Ad_wi;
					Addot_wi = this.body.Addot_wi;
					Addot_wp = this.parent.body.Addot_wi;
					Addot_ip = -Ad_iw*(Addot_wi*Ad_iw*Ad_wp - Addot_wp);
					colsA = jointA.idxR;
					J(rows,colsA) = Ad_ip*J(rowsP,colsA);
					Jdot(rows,colsA) = Ad_ip*Jdot(rowsP,colsA) + Addot_ip*J(rowsP,colsA);
					jointA = jointA.parent;
				end
				if ~isempty(this.next)
					[J,Jdot] = this.next.computeJacobian(J,Jdot);
				end
			end
		end
		
		%%
		function [y,z] = computeJacProd(this,x,y,z)
			% Computes y=J*x and z=Jdot*x
			nm = redmax.Scene.countM();
			if nargout == 1
				if nargin == 2
					y = zeros(nm,1);
				end
				rows = this.body.idxM;
				Ad_ij = this.body.Ad_ij;
				y(rows) = Ad_ij*this.S*x(this.idxR);
				if ~isempty(this.parent)
					rowsP = this.parent.body.idxM;
					Ad_ip = this.body.Ad_ip;
					y(rows) = y(rows) + Ad_ip*y(rowsP);
				end
				if ~isempty(this.next)
					y = this.next.computeJacProd(x,y);
				end
			else
				if nargin == 2
					y = zeros(nm,1);
					z = zeros(nm,1);
				end
				rows = this.body.idxM;
				Ad_ij = this.body.Ad_ij;
				y(rows) = Ad_ij*this.S*x(this.idxR);
				z(rows) = Ad_ij*this.Sdot*x(this.idxR);
				if ~isempty(this.parent)
					rowsP = this.parent.body.idxM;
					Ad_ip = this.body.Ad_ip;
					Ad_iw = this.body.Ad_iw;
					Ad_wp = this.parent.body.Ad_wi;
					Addot_wi = this.body.Addot_wi;
					Addot_wp = this.parent.body.Addot_wi;
					Addot_ip = -Ad_iw*(Addot_wi*Ad_iw*Ad_wp - Addot_wp);
					y(rows) = y(rows) + Ad_ip*y(rowsP);
					z(rows) = z(rows) + Ad_ip*z(rowsP) + Addot_ip*y(rowsP);
				end
				if ~isempty(this.next)
					[y,z] = this.next.computeJacProd(x,y,z);
				end
			end
		end
		
		%%
		function x = computeJacTransProd(this,y,x)
			% Computes x = J'*y
			nr = redmax.Scene.countR();
			if nargin == 2
				x = zeros(nr,1);
			end
			yi = y(this.body.idxM);
			for k = 1 : length(this.children)
				yi = yi + this.children{k}.alpha;
			end
			this.alpha = this.body.Ad_ip'*yi;
			x(this.idxR) = (this.body.Ad_ij*this.S)'*yi;
			if ~isempty(this.prev)
				x = this.prev.computeJacTransProd(y,x);
			end
		end
		
		%%
		function y = computeStiffnessProd(this,x,y)
			% Computes y = K*x
			nr = redmax.Scene.countR();
			if nargin == 2
				y = zeros(nr,1);
			end
			y(this.idxR) = -this.Kr*x(this.idxR); % negate to match rest of code
			if ~isempty(this.next)
				y = this.next.computeStiffnessProd(x,y);
			end
		end
		
		%%
		function y = computeDampingProd(this,x,y)
			% Computes y = D*x
			nr = redmax.Scene.countR();
			if nargin == 2
				y = zeros(nr,1);
			end
			y(this.idxR) = this.Dr*x(this.idxR);
			if ~isempty(this.next)
				y = this.next.computeDampingProd(x,y);
			end
		end
		
		%%
		function generateContacts(this)
			% Generates the local contact structures
			this.generateContacts_();
			for i = 1 : length(this.contacts)
				this.contacts{i}.a = 0;
			end
			if ~isempty(this.next)
				this.next.generateContacts();
			end
		end
		
		%%
		function fc = scatterContactForce(this,fc)
			% Computes the 3D world force that generates fc.

			% We know the overall maximal constraint force, fc, but we do
			% not know how to distribute these forces to each joint. For
			% example, if we have a double pendulum, fc is 12x1. The top 
			% 6x1 portion of this force is applied to one body and the
			% bottom 6x1 to the other body. For the parent body, the
			% constraint force is the sum of two joint reaction forces,
			% whereas for the child body, the constraint force is exactly
			% the joint reaction force. Assuming backward traversal, we can
			% process each joint/body one by one and subtracting the joint
			% reaction force from its parent before continuing the backward
			% traversal.
			this.body.fc_i = fc(this.body.idxM);
			%this.computeContactForce_(fci);
			% Subtract from this body's force
			fc(this.body.idxM) = zeros(6,1);
			% Subtract negated force from parent body's force
			if ~isempty(this.parent)
				Ad = se3.Ad(this.body.E_iw * this.parent.body.E_wi);
				fc(this.prev.body.idxM) = fc(this.prev.body.idxM) + Ad'*this.body.fc_i;
			end
			% Backward traversal
			if ~isempty(this.prev)
				fc = this.prev.scatterContactForce(fc);
			end
		end
		
		%%
		function computeContactMultiplier(this,h,SPreg)
			% Computes the contact Lagrange multiplier
			nc = length(this.contacts);
			if isempty(this.parent)
				Minv = diag(1./this.body.I_i);
				fc = this.body.fc_i;
				N = zeros(nc,6);
			else
				E_ip = this.body.E_iw*this.parent.body.E_wi;
				Minv = diag(1./[this.body.I_i;this.parent.body.I_i]);
				fc(1:6,1) = this.body.fc_i;
				fc(7:12,1) = -se3.Ad(E_ip)'*this.body.fc_i;
				N = zeros(nc,12);
			end
			for i = 1 : nc
				contact = this.contacts{i};
				nor_i = contact.nor_i;
				pos_i = contact.pos_i;
				N(i,1:6) = nor_i'*se3.Gamma(pos_i);
				if ~isempty(this.parent)
					E_pi = this.parent.body.E_iw*this.body.E_wi;
					nor_p = E_pi(1:3,:)*[nor_i;0];
					pos_p = E_pi(1:3,:)*[pos_i;1];
					N(i,7:12) = -nor_p'*se3.Gamma(pos_p);
				end
			end
			% Since fa = -N'*a/h, we need to multiply by h to get a
			H = N*Minv*N' + SPreg*eye(nc);
			H = 0.5*(H + H');
			f = (N*Minv*fc)*h;
			a = H\f;
			for i = 1 : nc
				this.contacts{i}.a = a(i);
			end
			
			if ~isempty(this.next)
				this.next.computeContactMultiplier(h,SPreg);
			end
		end
		
		%%
		function T = computeTangentMatrix(this,T)
			if nargin == 1
				nt = redmax.Scene.countT();
				nm = redmax.Scene.countM();
				T = zeros(nt,nm);
			end
			T = this.computeTangentMatrix_(T);
			if ~isempty(this.next)
				T = this.next.computeTangentMatrix(T);
			end
		end
		
		%%
		function [bl,bu,idx] = computeFrictionLimits(this,mu,SPathresh,bl,bu,idx)
			if nargin == 3
				nt = redmax.Scene.countT();
				bl = zeros(nt,1);
				bu = zeros(nt,1);
				idx = [];
			end
			[bl,bu,idx] = this.computeFrictionLimits_(mu,SPathresh,bl,bu,idx);
			if ~isempty(this.next)
				[bl,bu,idx] = this.next.computeFrictionLimits(mu,SPathresh,bl,bu,idx);
			end
		end
		
		%%
		function [T,V] = computeEnergies(this,grav,T,V)
			% Computes kinetic and potential energies
			if nargin == 2
				T = 0;
				V = 0;
			end
			[T,V] = this.body.computeEnergies(grav,T,V);
			dq = this.q - this.qInit;
			V = V + 0.5*this.Kr*(dq'*dq);
			if ~isempty(this.next)
				[T,V] = this.next.computeEnergies(grav,T,V);
			end
		end
		
		%%
		function draw(this)
			%se3.drawAxis(this.E_wj,this.body.getAxisSize());
			this.draw_();
			if ~isempty(this.next)
				this.next.draw();
			end
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function reparam_(this) %#ok<MANU>
		end
		
		%%
		function update_(this) %#ok<MANU>
		end
		
		%%
		function generateContacts_(this) %#ok<MANU>
		end
		
		%%
		function T = computeTangentMatrix_(this,T) %#ok<INUSL>
		end
		
		%%
		function [bl,bu,idx] = computeFrictionLimits_(this,mu,SPathresh,bl,bu,idx) %#ok<INUSL>
		end
		
		%%
		function draw_(this) %#ok<MANU>
		end
	end
end
