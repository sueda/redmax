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
		Kr       % Joint stiffness
		Dr       % Joint damping
		S        % Jacobian
		Sdot     % dS/dt
		I_j      % Inertia at the joint
		V        % Twist at parent joint
		Vdot     % Acceleration at parent joint
		E_pj     % Transform of this joint wrt parent joint
		E_pj0    % Transform when q is zero (initial transform)
		E_jp     % Transform of parent joint wrt this joint
		Ad_jp    % Adjoint of E_jp
		E_wj     % Transform of this joint wrt world
		next     % Forward recursive ordering
		prev     % Reverse recursive ordering
		idxR     % Reduced indices
	end
	
	%%
	properties (Access = protected)
		Q % Transformation matrix applied about the joint
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
			this.Kr = 0;
			this.Dr = 0;
			this.S = zeros(6,ndof);
			this.Sdot = zeros(6,ndof);
			this.I_j = eye(6);
			this.V = zeros(6,1);
			this.Vdot = zeros(6,1);
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
			% Gets traversal order. Given the tree structure (parent,
			% children) of the joints, this function returns an array
			% representation so that the parent comes before all of its
			% children.
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
			% which can also be computed as:
			%   this.I_j = this.body.Ad_ij'*diag(this.body.I_i)*this.body.Ad_ij;
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
			% Update attached body
			this.body.update();
			if ~isempty(this.next)
				this.next.update();
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
			rows = this.idxR;
			% We'll add the joint torque here rather than having a
			% separate function.
			fr(rows) = fr(rows) + this.tau - this.Kr*(this.q - this.qInit);
			Kr(rows,rows) = Kr(rows,rows) - this.Kr*eye(this.ndof);
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
			rows = this.idxR;
			fr(rows) = fr(rows) - this.Dr*this.qdot;
			Dr(rows,rows) = Dr(rows,rows) + this.Dr*eye(this.ndof);
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
				% Only compute J
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
				% Compute J and Jdot
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
		function update_(this) %#ok<MANU>
		end
		
		%%
		function draw_(this) %#ok<MANU>
		end
	end
end
