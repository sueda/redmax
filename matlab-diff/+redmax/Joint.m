classdef (Abstract) Joint < handle
	%Joint A generic joint between two bodies
	%   A joint is defined between a parent and the child. The DOF, q, of
	%   the joint is the relative displacement of the child wrt the parent.
	
	%%
	properties
		name      % Optional name
		body      % Attached body
		parent    % Parent joint
		children  % Children joints
		ndof      % Number of DOF
		q         % Position
		qdot      % Velocity
		q0        % Last position
		qdot0     % Last velocity
		q1        % Last position (for BDF2, q1 is k, and q0 is k-1)
		qdot1     % Last velocity (for BDF2, qdot1 is k, and qdot0 is k-1)
		qRest     % Rest position (for joint stiffness)
		tau       % Joint torque
		stiffness % Joint stiffness
		damping   % Joint damping
		I_j       % Inertia at the joint
		V         % Twist at parent joint
		E_pj      % Transform of this joint wrt parent joint
		E_jp      % Transform of parent joint wrt this joint
		E0_pj     % Transform of this joint wrt parent joint when q=0
		E0_jp     % Transform of parent joint wrt this joint when q=0
		A_jp      % Adjoint of E_jp
		E_wj      % Transform of this joint wrt world
		next      % Forward recursive ordering
		prev      % Reverse recursive ordering
		idxR      % Reduced indices
		
		% Quantities to be computed by subclasses in update()
		Q         % Transformation matrix applied about the joint
		A         % Adjoint of Q
		dAdq      % Derivative of A wrt q
		Adot      % Derivative of A wrt t
		dAdotdq   % Derivative of Adot wrt q
		S         % Jacobian
		dSdq      % Derivative of S wrt q
		Sdot      % Derivative of S wrt t
		dSdotdq   % Derivative of Sdot wrt q
		invQ      % inv(Q)
		invA      % inv(A)
	end
	
	%%
	methods
		%%
		function this = Joint(parent,body,ndof)
			pname = 'NULL';
			bname = 'NULL';
			if ~isempty(parent)
				if ~isempty(parent.body)
					pname = parent.body.name;
				end
			end
			if ~isempty(body)
				bname = body.name;
			end
			this.name = [pname,'-',bname];
			this.parent = parent;
			this.body = body;
			this.ndof = ndof;
			this.q = zeros(ndof,1);
			this.qdot = zeros(ndof,1);
			this.q0 = zeros(ndof,1);
			this.qdot0 = zeros(ndof,1);
			this.q1 = zeros(ndof,1);
			this.qdot1 = zeros(ndof,1);
			this.tau = zeros(ndof,1);
			this.stiffness = 0;
			this.damping = 0;
			this.S = zeros(6,ndof);
			this.Sdot = zeros(6,ndof);
			this.I_j = eye(6);
			this.V = zeros(6,1);
			body.joint = this;
			if ~isempty(parent)
				parent.children{end+1} = this;
			end
		end
		
		%%
		function setJointTransform(this,E)
			% Sets the transform of this joint wrt parent joint
			this.E0_pj = E;
			this.E0_jp = se3.inv(E);
		end
		
		%%
		function setStiffness(this,stiffness)
			% Sets this joint's linear stiffness
			this.stiffness = stiffness;
		end
		
		%%
		function setDamping(this,damping)
			% Sets this joint's linear velocity damping
			this.damping = damping;
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
			this.qRest = this.q;
		end
		
		%%
		function computeInertia(this)
			m = this.body.I_i(4);
			R = this.body.E0_ji(1:3,1:3);
			p = this.body.E0_ji(1:3,4);
			pBrac = se3.brac(p);
			Ic = diag(this.body.I_i(1:3));
			this.I_j = [R*Ic*R'+m*(pBrac'*pBrac),m*pBrac;m*pBrac',m*eye(3)];
			% which can also be computed as:
			%   this.I_j = this.body.A_ij'*diag(this.body.I_i)*this.body.A_ij;
		end
		
		%%
		function [qAll,qdotAll] = getQ(this,qAll,qdotAll)
			% Gathers q (and optionally qdot) into (a) global vector(s)
			nr = redmax.Scene.countR();
			if nargin == 1
				qAll = zeros(nr,1);
				if nargout == 2
					qdotAll = zeros(nr,1);
				end
			end
			qAll(this.idxR) = this.q;
			if nargout == 2
				qdotAll(this.idxR) = this.qdot;
			end
			% Go to the next joint
			if ~isempty(this.next)
				if nargout == 1
					qAll = this.next.getQ(qAll);
				else
					[qAll,qdotAll] = this.next.getQ(qAll,qdotAll);
				end
			end
		end
		
		%%
		function qdotAll = getQdot(this,qdotAll)
			% Gathers qdot into a global vector
			nr = redmax.Scene.countR();
			if nargin == 1
				qdotAll = zeros(nr,1);
			end
			qdotAll(this.idxR) = this.qdot;
			% Go to the next joint
			if ~isempty(this.next)
				qdotAll = this.next.getQdot(qdotAll);
			end
		end
		
		%%
		function setQ(this,qAll,qdotAll)
			% Scatters q (and optionally qdot) from (a) global vector(s)
			this.q(1:this.ndof) = qAll(this.idxR);
			if nargin == 3
				this.qdot(1:this.ndof) = qdotAll(this.idxR);
			end
			% Go to the next joint
			if ~isempty(this.next)
				if nargin == 2
					this.next.setQ(qAll);
				else
					this.next.setQ(qAll,qdotAll);
				end
			end
		end
		
		%%
		function setQdot(this,qdotAll)
			% Scatters qdot from a global vector
			this.qdot(1:this.ndof) = qdotAll(this.idxR);
			% Go to the next joint
			if ~isempty(this.next)
				this.next.setQdot(qdotAll);
			end
		end
		
		%%
		function [q0All,qdot0All] = getQ0(this,q0All,qdot0All)
			% Gathers q0 (and optionally qdot0) into (a) global vector(s)
			nr = redmax.Scene.countR();
			if nargin == 1
				q0All = zeros(nr,1);
				if nargout == 2
					qdot0All = zeros(nr,1);
				end
			end
			q0All(this.idxR) = this.q0;
			if nargout == 2
				qdot0All(this.idxR) = this.qdot0;
			end
			% Go to the next joint
			if ~isempty(this.next)
				if nargout == 1
					q0All = this.next.getQ0(q0All);
				else
					[q0All,qdot0All] = this.next.getQ0(q0All,qdot0All);
				end
			end
		end
		
		%%
		function qdot0All = getQdot0(this,qdot0All)
			% Gathers qdot0 into a global vector
			nr = redmax.Scene.countR();
			if nargin == 1
				qdot0All = zeros(nr,1);
			end
			qdot0All(this.idxR) = this.qdot0;
			% Go to the next joint
			if ~isempty(this.next)
				qdot0All = this.next.getQdot0(qdot0All);
			end
		end
		
		%%
		function setQ0(this,q0All,qdot0All)
			% Scatters q0 (and optionally qdot0) from (a) global vector(s)
			this.q0(1:this.ndof) = q0All(this.idxR);
			if nargin == 3
				this.qdot0(1:this.ndof) = qdot0All(this.idxR);
			end
			% Subclasses may need to save other things
			this.setAux0_();
			% Go to the next joint
			if ~isempty(this.next)
				if nargin == 2
					this.next.setQ0(q0All);
				else
					this.next.setQ0(q0All,qdot0All);
				end
			end
		end
		
		%%
		function setQdot0(this,qdot0All)
			% Scatters qdot0 from a global vector
			this.qdot0(1:this.ndof) = qdot0All(this.idxR);
			% Go to the next joint
			if ~isempty(this.next)
				this.next.setQdot0(qdot0All);
			end
		end
		
		%%
		function [q1All,qdot1All] = getQ1(this,q1All,qdot1All)
			% Gathers q1 (and optionally qdot1) into (a) global vector(s)
			nr = redmax.Scene.countR();
			if nargin == 1
				q1All = zeros(nr,1);
				if nargout == 2
					qdot1All = zeros(nr,1);
				end
			end
			q1All(this.idxR) = this.q1;
			if nargout == 2
				qdot1All(this.idxR) = this.qdot1;
			end
			% Go to the next joint
			if ~isempty(this.next)
				if nargout == 1
					q1All = this.next.getQ1(q1All);
				else
					[q1All,qdot1All] = this.next.getQ1(q1All,qdot1All);
				end
			end
		end
		
		%%
		function qdot1All = getQdot1(this,qdot1All)
			% Gathers qdot1 into a global vector
			nr = redmax.Scene.countR();
			if nargin == 1
				qdot1All = zeros(nr,1);
			end
			qdot1All(this.idxR) = this.qdot1;
			% Go to the next joint
			if ~isempty(this.next)
				qdot1All = this.next.getQdot1(qdot1All);
			end
		end
		
		%%
		function setQ1(this,q1All,qdot1All)
			% Scatters q1 (and optionally qdot1) from (a) global vector(s)
			this.q1(1:this.ndof) = q1All(this.idxR);
			if nargin == 3
				this.qdot1(1:this.ndof) = qdot1All(this.idxR);
			end
			% Subclasses may need to save other things
			this.setAux1_();
			% Go to the next joint
			if ~isempty(this.next)
				if nargin == 2
					this.next.setQ1(q1All);
				else
					this.next.setQ1(q1All,qdot1All);
				end
			end
		end
		
		%%
		function setQdot1(this,qdot1All)
			% Scatters qdot1 from a global vector
			this.qdot1(1:this.ndof) = qdot1All(this.idxR);
			% Go to the next joint
			if ~isempty(this.next)
				this.next.setQdot1(qdot1All);
			end
		end
		
		%%
		function reparam(this)
			% Subclasses may need to reparameterize q and qdot
			this.reparam_();
			% Go to the next joint
			if ~isempty(this.next)
				this.next.reparam();
			end
		end
		
		%%
		function update(this,deriv)
			if nargin == 1
				deriv = true;
			end
			% Reset
			n = this.ndof;
			this.Q = eye(4,4);
			this.A = eye(6,6);
			this.Adot = zeros(6,6);
			this.S = zeros(6,n);
			this.Sdot = zeros(6,n);
			if deriv
				this.dAdq = zeros(6,6,n);
				this.dAdotdq = zeros(6,6,n);
				this.dSdq = zeros(6,n,n);
				this.dSdotdq = zeros(6,n,n);
			end
			% Updates this joint and the attached body
			this.update_(deriv);
			this.invQ = se3.inv(this.Q);
			this.invA = se3.Ad(this.invQ);
			% Transforms and adjoints
			if isempty(this.E0_pj)
				this.E_pj = this.Q;
			else
				this.E_pj = this.E0_pj*this.Q;
			end
			this.E_jp = se3.inv(this.E_pj);
			this.A_jp = se3.Ad(this.E_jp);
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
				this.V = this.V + this.A_jp*this.parent.V;
			end
			% Update attached body
			if ~isempty(this.body)
				this.body.update();
			end
			% Go to the next joint
			if ~isempty(this.next)
				this.next.update();
			end
		end
		
		%%
		function [fr,Kr,Dr] = computeForce(this,fr,Kr,Dr)
			% Computes reduced force vector and matrix
			nr = redmax.Scene.countR();
			if nargout == 1
				% Just the force
				if nargin == 1
					fr = zeros(nr,1);
				end
				fr = this.computeForce_(fr);
				rows = this.idxR;
				% Joint torque
				fr(rows) = fr(rows) + this.tau + this.stiffness*(this.qRest - this.q) - this.damping*this.qdot;
				% Go to the next joint
				if ~isempty(this.next)
					fr = this.next.computeForce(fr);
				end
			else
				% Both force and derivatives
				if nargin == 1
					fr = zeros(nr,1);
					Kr = zeros(nr);
					Dr = zeros(nr);
				end
				[fr,Kr,Dr] = this.computeForce_(fr,Kr,Dr);
				rows = this.idxR;
				% Joint torque
				fr(rows) = fr(rows) + this.tau + this.stiffness*(this.qRest - this.q) - this.damping*this.qdot;
				I = eye(this.ndof);
				Kr(rows,rows) = Kr(rows,rows) - this.stiffness*I;
				Dr(rows,rows) = Dr(rows,rows) - this.damping*I;
				% Go to the next joint
				if ~isempty(this.next)
					[fr,Kr,Dr] = this.next.computeForce(fr,Kr,Dr);
				end
			end
		end
		
		%%
		function [J,Jdot,dJdq,dJdotdq] = computeJacobian(this,J,Jdot,dJdq,dJdotdq)
			% Computes the Jacobian and its derivatives
			nm = redmax.Scene.countM();
			nr = redmax.Scene.countR();
			if nargout == 2
				% Only J and Jdot (takes O(n^2) time)
				if nargin == 1
					J = zeros(nm,nr);
					Jdot = zeros(nm,nr);
				end
				invQ = this.invQ; %#ok<*PROPLC>
				Adot = this.Adot;
				S = this.S;
				Sdot = this.Sdot;
				idxmI = this.body.idxM;
				idxrI = this.idxR;
				E0_BiJi = this.body.E0_ij; % this body - this joint
				A0_BiJi = this.body.A0_ij;
				J   (idxmI,idxrI) = A0_BiJi*S;
				Jdot(idxmI,idxrI) = A0_BiJi*Sdot;
				if ~isempty(this.parent)
					idxmP = this.parent.body.idxM;
					E0_JpBp = this.parent.body.E0_ji; % parent joint - parent body
					E0_JiJp = this.E0_jp;             % this joint   - parent joint
					E0_JiBp = E0_JiJp*E0_JpBp;        % this joint   - parent body
					E_BiBp = E0_BiJi*invQ*E0_JiBp;    % this body    - parent body
					A_BiBp = se3.Ad(E_BiBp);
					Aleft  = -se3.Ad(E0_BiJi*invQ);
					Aright =  se3.Ad(invQ*E0_JiBp);
					Adot_BiBp = Aleft*Adot*Aright;
					jointA = this.parent;
					while ~isempty(jointA)
						idxrA = jointA.idxR;
						JPA    = J   (idxmP,idxrA);
						JdotPA = Jdot(idxmP,idxrA);
						J   (idxmI,idxrA) = A_BiBp*JPA;
						Jdot(idxmI,idxrA) = A_BiBp*JdotPA + Adot_BiBp*JPA;
						jointA = jointA.parent;
					end
				end
				% Go to the next joint
				if ~isempty(this.next)
					[J,Jdot] = this.next.computeJacobian(J,Jdot);
				end
			else
				% J, Jdot, dJdq, dJdotdq  (takes O(n^3) time)
				if nargin == 1
					J = zeros(nm,nr);
					Jdot = zeros(nm,nr);
					dJdq = zeros(nm,nr,nr);
					dJdotdq = zeros(nm,nr,nr);
				end
				invQ = this.invQ; %#ok<*PROPLC>
				invA = this.invA;
				dAdq = this.dAdq;
				Adot = this.Adot;
				dAdotdq = this.dAdotdq;
				S = this.S;
				dSdq = this.dSdq;
				Sdot = this.Sdot;
				dSdotdq = this.dSdotdq;
				idxmI = this.body.idxM;
				idxrI = this.idxR;
				E0_BiJi = this.body.E0_ij; % this body - this joint
				A0_BiJi = this.body.A0_ij;
				J   (idxmI,idxrI) = A0_BiJi*S;
				Jdot(idxmI,idxrI) = A0_BiJi*Sdot;
				for ii = 1 : this.ndof
					dJdq   (idxmI,idxrI,idxrI(ii)) = A0_BiJi*   dSdq(:,:,ii);
					dJdotdq(idxmI,idxrI,idxrI(ii)) = A0_BiJi*dSdotdq(:,:,ii);
				end
				if ~isempty(this.parent)
					idxmP = this.parent.body.idxM;
					E0_JpBp = this.parent.body.E0_ji; % parent joint - parent body
					E0_JiJp = this.E0_jp;             % this joint   - parent joint
					E0_JiBp = E0_JiJp*E0_JpBp;        % this joint   - parent body
					E_BiBp = E0_BiJi*invQ*E0_JiBp;    % this body    - parent body
					A_BiBp = se3.Ad(E_BiBp);
					dAdq_BiBp    = zeros(6,6,this.ndof);
					dAdotdq_BiBp = zeros(6,6,this.ndof);
					Aleft  = -se3.Ad(E0_BiJi*invQ);
					Aright =  se3.Ad(invQ*E0_JiBp);
					Adot_BiBp = Aleft*Adot*Aright;
					for ii = 1 : this.ndof
						dAdq_ii    =    dAdq(:,:,ii);
						dAdotdq_ii = dAdotdq(:,:,ii);
						tmp1 = dAdq_ii*invA*Adot;
						tmp2 = Adot*invA*dAdq_ii;
						dAdq_BiBp   (:,:,ii) = Aleft*dAdq_ii*Aright;
						dAdotdq_BiBp(:,:,ii) = Aleft*(dAdotdq_ii - tmp1 - tmp2)*Aright;
					end
					jointA = this.parent;
					while ~isempty(jointA)
						idxrA = jointA.idxR;
						JPA    = J   (idxmP,idxrA);
						JdotPA = Jdot(idxmP,idxrA);
						J   (idxmI,idxrA) = A_BiBp*JPA;
						Jdot(idxmI,idxrA) = A_BiBp*JdotPA + Adot_BiBp*JPA;
						for ii = 1 : length(idxrI)
							dAdq_BiBp_ii    = dAdq_BiBp   (:,:,ii);
							dAdotdq_BiBp_ii = dAdotdq_BiBp(:,:,ii);
							dJdq   (idxmI,idxrA,idxrI(ii)) = dAdq_BiBp_ii*JPA;
							dJdotdq(idxmI,idxrA,idxrI(ii)) = dAdq_BiBp_ii*JdotPA + dAdotdq_BiBp_ii*JPA;
						end
						jointK = this.parent;
						while ~isempty(jointK)
							idxrK = jointK.idxR;
							for kk = 1 : length(idxrK)
								dJdqPAK    = dJdq   (idxmP,idxrA,idxrK(kk));
								dJdotdqPAK = dJdotdq(idxmP,idxrA,idxrK(kk));
								dJdq   (idxmI,idxrA,idxrK(kk)) = A_BiBp*dJdqPAK;
								dJdotdq(idxmI,idxrA,idxrK(kk)) = A_BiBp*dJdotdqPAK + Adot_BiBp*dJdqPAK;
							end
							jointK = jointK.parent;
						end
						jointA = jointA.parent;
					end
				end
				% Go to the next joint
				if ~isempty(this.next)
					[J,Jdot,dJdq,dJdotdq] = this.next.computeJacobian(J,Jdot,dJdq,dJdotdq);
				end
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
			dq = this.q - this.qRest;
			V = V + 0.5*this.stiffness*(dq'*dq);
			% Go to the next joint
			if ~isempty(this.next)
				[T,V] = this.next.computeEnergies(grav,T,V);
			end
		end
		
		%%
		function draw(this)
			%se3.drawAxis(this.E_wj,this.body.getAxisSize());
			this.draw_();
			% Go to the next joint
			if ~isempty(this.next)
				this.next.draw();
			end
		end

		%%
		function test(this)
			this.update();
			n = this.ndof;
			q = this.q; %#ok<*PROP>
			Q = this.Q;
			A = this.A;
			dAdq = this.dAdq;
			Adot = this.Adot;
			dAdotdq = this.dAdotdq;
			S = this.S;
			dSdq = this.dSdq;
			Sdot = this.Sdot;
			dSdotdq = this.dSdotdq;
			
			fprintf('=== %s ===\n',class(this));
			h = sqrt(eps);

			% Check dAdq, dAdotdq, dSdq, Sdotdq, and S
			dAdq_ = zeros(6,6,n);
			dSdq_ = zeros(6,n,n);
			dAdotdq_ = zeros(6,6,n);
			dSdotdq_ = zeros(6,n,n);
			S_ = zeros(6,n);
			for j = 1 : n
				q_ = q;
				q_(j) = q_(j) + h;
				this.q = q_;
				this.update();
				A_ = this.A;
				S_ = this.S;
				Q_ = this.Q;
				Adot_ = this.Adot;
				Sdot_ = this.Sdot;
				this.q = q;
				this.update();
				dAdq_(:,:,j) = (A_ - A)/h;
				dSdq_(:,:,j) = (S_ - S)/h;
				dQdq_ = (Q_ - Q)/h;
				dAdotdq_(:,:,j) = (Adot_ - Adot)/h;
				dSdotdq_(:,:,j) = (Sdot_ - Sdot)/h;
				S_(:,j) = se3.unbrac(Q\dQdq_);
			end
			redmax.Scene.printError('dA/dq',dAdq,dAdq_);
			redmax.Scene.printError('dS/dq',dSdq,dSdq_);
			redmax.Scene.printError('dAdot/dq',dAdotdq,dAdotdq_);
			redmax.Scene.printError('dSdot/dq',dSdotdq,dSdotdq_);
			redmax.Scene.printError('S',S,S_);

			% Check Adot and Sdot
			Adot = zeros(6);
			Sdot_ = zeros(6,n); % another way to compute Sdot
			for j = 1 : n
				Adot = Adot + dAdq(:,:,j)*this.qdot(j);
				Sdot_ = Sdot_ + dSdq(:,:,j)*this.qdot(j);
			end
			redmax.Scene.printError('Sdot',Sdot,Sdot_);
			q = this.q;
			q_ = q + h*this.qdot;
			this.q = q_;
			this.update();
			A_ = this.A;
			S_ = this.S;
			this.q = q;
			this.update();
			Sdot_ = (S_ - S)/h;
			Adot_ = (A_ - A)/h;
			redmax.Scene.printError('Sdot',Sdot,Sdot_);
			redmax.Scene.printError('Adot',Adot,Adot_);
			
			% Check Adot_BiBp, dAdq_BiBp, and dAdotdq_BiBp
			if ~isempty(this.parent)
				invQ = se3.inv(Q);
				E_JiBi = this.body.E0_ji;         % this joint   - this body
				E_BiJi = se3.inv(E_JiBi);         % this body    - this joint
				E_JpBp = this.parent.body.E0_ji;  % parent joint - parent body
				E_JiJp = se3.inv(this.E_pj);      % this joint   - parent joint
				E_JiBp = E_JiJp*E_JpBp;           % this joint   - parent body
				E_BiBp = E_BiJi*invQ*E_JiBp;      % this body    - parent body
				A_BiBp = se3.Ad(E_BiBp);
				dAdq_BiBp    = zeros(6,6,n);
				dAdotdq_BiBp = zeros(6,6,n);
				Aleft  = -se3.Ad(E_BiJi*invQ);
				Aright =  se3.Ad(invQ*E_JiBp);
				invA = se3.Ad(invQ);
				Adot_BiBp = Aleft*Adot*Aright;
				for ii = 1 : n
					dAdqi    =    dAdq(:,:,ii);
					dAdotdqi = dAdotdq(:,:,ii);
					tmp1 = dAdqi*invA*Adot;
					tmp2 = Adot*invA*dAdqi;
					dAdq_BiBp   (:,:,ii) = Aleft*dAdqi*Aright;
					dAdotdq_BiBp(:,:,ii) = Aleft*(dAdotdqi - tmp1 - tmp2)*Aright;
				end
				% Check...
				q = this.q;
				q_ = q + h*this.qdot;
				this.q = q_;
				this.update();
				Q_ = this.Q;
				this.q = q;
				this.update();
				invQ_ = se3.inv(Q_);
				E_BiBp_ = E_BiJi*invQ_*E_JiBp;
				A_BiBp_ = se3.Ad(E_BiBp_);
				Adot_BiBp_ = (A_BiBp_ - A_BiBp)/h;
				redmax.Scene.printError('Adot_BiBp',Adot_BiBp_,Adot_BiBp);
				dAdq_BiBp_    = zeros(6,6,n);
				dAdotdq_BiBp_ = zeros(6,6,n);
				for j = 1 : n
					q = this.q;
					q_ = q;
					q_(j) = q_(j) + h;
					this.q = q_;
					this.update();
					Q_ = this.Q;
					Adot_ = this.Adot;
					this.q = q;
					this.update();
					invQ_ = se3.inv(Q_);
					E_BiBp_ = E_BiJi*invQ_*E_JiBp;
					A_BiBp_ = se3.Ad(E_BiBp_);
					Aleft_  = -se3.Ad(E_BiJi*invQ_);
					Aright_ =  se3.Ad(invQ_*E_JiBp);
					Adot_BiBp_ = Aleft_*Adot_*Aright_;
					dAdq_BiBp_(:,:,j) = (A_BiBp_ - A_BiBp)/h;
					dAdotdq_BiBp_(:,:,j) = (Adot_BiBp_ - Adot_BiBp)/h;
				end
				redmax.Scene.printError('dAdq_BiBp_',dAdq_BiBp,dAdq_BiBp_);
				redmax.Scene.printError('dAdotdq_BiBp_',dAdotdq_BiBp,dAdotdq_BiBp_);
			end
			% Check inverse Euler
			if isa(this,'redmax.JointSpherical')
				this.testEuler(); %#ok<MCNPN>
			end
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function reparam_(this) %#ok<MANU>
			% Subclasses may need to reparameterize q and qdot
		end
		
		%%
		function setAux0_(this) %#ok<MANU>
			% Subclasses may need to save other things
		end
		
		%%
		function setAux1_(this) %#ok<MANU>
			% Subclasses may need to save other things
		end
		
		%%
		function update_(this,deriv) %#ok<INUSD>
			% Compute Q,A,Adot,S,Sdot
			% if deriv, also compute dAdq,dAdotdq,dSdq,dSdotdq
		end
		
		%%
		function [fr,Kr,Dr] = computeForce_(this,fr,Kr,Dr) %#ok<INUSL>
			% Subclass specific force
		end
		
		%%
		function draw_(this) %#ok<MANU>
			% Subclass specific draw
		end
	end
end
