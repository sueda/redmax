classdef (Abstract) Body < handle
	%Body A rigid body connected through a joint to a parent
	
	%%
	properties
		name     % Optional name
		color    % Color for rendering
		density  % Mass/volume
		I_i      % Inertia at body center
		E_ji     % Where the body is wrt joint
		E_ij     % Where the joint is wrt body
		E_wi     % Where the body is wrt world
		E_iw     % Where the world is wrt body
		E_ip     % Where the parent is wrt body
		Ad_ji    % Adjoint of E_ji
		Ad_ij    % Adjoint of E_ij
		Ad_iw    % Adjoint of E_iw
		Ad_wi    % Adjoint of E_wi
		Ad_ip    % Adjoint of E_ip
		Addot_wi % Adjoint dot of E_wi, used by the Jacobian
		phi      % Twist at body center
		phidot   % Acceleration at body center
		damping  % Viscous damping
		joint    % Joint to parent
		idxM     % Maximal indices
		next     % Next body in traversal order
	end
	
	%%
	methods
		%%
		function this = Body(density)
			global countB countCM CM;
			this.name = ['body',num2str(countB)];
			this.color = CM(mod(countCM-1,size(CM,1))+1,:);
			this.density = density;
			this.I_i = eye(6);
			this.E_ji = eye(4);
			this.E_ij = eye(4);
			this.E_wi = eye(4);
			this.E_iw = eye(4);
			this.E_ip = eye(4);
			this.Ad_ji = eye(6);
			this.Ad_ij = eye(6);
			this.Ad_iw = eye(6);
			this.Ad_wi = eye(6);
			this.Ad_ip = eye(6);
			this.Addot_wi = eye(6);
			this.phi = zeros(6,1);
			this.phidot = zeros(6,1);
			this.damping = 0;
			this.joint = {};
			countB = countB + 1;
			countCM = countCM + 1;
		end
		
		%%
		function setBodyTransform(this,E)
			% Sets the transform of this body wrt parent joint
			this.E_ji = E;
			this.E_ij = se3.inv(this.E_ji);
			this.Ad_ji = se3.Ad(this.E_ji);
			this.Ad_ij = se3.Ad(this.E_ij);
		end
		
		%%
		function setDamping(this,damping)
			% Sets viscous damping
			this.damping = damping;
		end
		
		%%
		function countDofs(this)
			% Counts maximal DOFs
			nm = redmax.Scene.countM();
			this.idxM = nm + (1:6);
			nm = nm + 6;
			redmax.Scene.countM(nm);
		end
		
		%%
		function computeInertia(this)
			% Computes inertia at body and joint
			this.computeInertia_();
			this.joint.computeInertia();
		end
		
		%%
		function update(this)
			% Updates this body's transforms and velocities
			this.E_wi = this.joint.E_wj*this.E_ji;
			this.E_iw = se3.inv(this.E_wi);
			this.Ad_wi = se3.Ad(this.E_wi);
			this.Ad_iw = se3.Ad(this.E_iw);
			this.E_ip = eye(4);
			if ~isempty(this.joint.parent)
				this.E_ip = this.E_iw*this.joint.parent.body.E_wi;
			end
			this.Ad_ip = se3.Ad(this.E_ip);
			% Body velocity
			this.phi = this.Ad_ij*this.joint.V;
			this.Addot_wi = se3.Addot(this.E_wi,this.phi);
			this.phidot = this.Ad_ij*this.joint.Vdot;
		end
		
		%%
		function [M,f] = computeMassGrav(this,grav,M,f)
			% Computes maximal mass matrix and force vector
			if nargin == 2
				nm = redmax.Scene.countM();
				M = zeros(nm,nm);
				f = zeros(nm,1);
			end
			rows = this.idxM;
			M_i = diag(this.I_i);
			M(rows,rows) = M_i;
			fcor = se3.ad(this.phi)'*M_i*this.phi;
			R_wi = this.E_wi(1:3,1:3);
			R_iw = R_wi';
			fgrav = zeros(6,1);
			fgrav(4:6) = M_i(4,4)*R_iw*grav; % wrench in body space
			f(rows) = fcor + fgrav;
			% Joint torque
			% This is how we would apply a joint torque using maximal
			% coordinates. It's much easier in reduced, so we'll do that
			% instad. See Joint.computeForce().
			%
			%	tau = this.joint.S*(this.joint.tau - this.joint.K*this.joint.q);
			%	f(rows) = f(rows) + this.Ad_ji'*tau;
			%	% Also apply to parent
			%	if ~isempty(this.joint.parent)
			%		parent = this.joint.parent.body;
			%		rowsP = parent.idxM;
			%		E_jp = this.E_ji*this.E_iw*parent.E_wi; % this joint -> parent body
			%		f(rowsP) = f(rowsP) - se3.Ad(E_jp)'*tau;
			%	end
			%
			if ~isempty(this.next)
				[M,f] = this.next.computeMassGrav(grav,M,f);
			end
		end
		
		%%
		function [f,D] = computeForceDamping(this,f,D)
			% Computes maximal damping force vector and matrix
			nm = redmax.Scene.countM();
			if nargin == 1
				f = zeros(nm,1);
				D = zeros(nm);
			elseif nargin == 2
				D = zeros(nm);
			end
			if this.damping > 0
				rows = this.idxM;
				fi = -this.damping*this.phi;
				Di = this.damping*eye(6);
				f(rows) = f(rows) + fi;
				D(rows,rows) = D(rows,rows) + Di;
			end
			if ~isempty(this.next)
				[f,D] = this.next.computeForceDamping(f,D);
			end
		end
		
		%%
		function draw(this)
			%s = this.getAxisSize();
			%if s > 0
			%	se3.drawAxis(this.E_wi,s);
			%end
			%text(this.E_wi(1,4),this.E_wi(2,4),this.E_wi(3,4),this.name);
			this.draw_();
			if ~isempty(this.next)
				this.next.draw();
			end
		end
		
		%%
		function s = getAxisSize(this) %#ok<MANU>
			s = 1;
		end
	end
	
	%%
	methods (Abstract)
		%%
		draw_(this)
		computeInertia_(this)
	end
end
