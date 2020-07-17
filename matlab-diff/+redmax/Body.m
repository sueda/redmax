classdef (Abstract) Body < handle
	%Body A rigid body connected through a joint to a parent
	
	%%
	properties
		name     % Optional name
		color    % Color for rendering
		density  % Mass/volume
		damping  % Viscous damping
		I_i      % Inertia at body center
		E0_ji    % Where the body is wrt joint (fixed)
		E0_ij    % Where the joint is wrt body (fixed)
		E_wi     % Where the body is wrt world
		E_iw     % Where the world is wrt body
		E_ip     % Where the parent is wrt body
		A0_ij    % Adjoint of E0_ij (fixed)
		phi      % Twist at body center
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
			this.damping = 0;
			this.I_i = eye(6);
			this.E0_ji = eye(4);
			this.E0_ij = eye(4);
			this.E_wi = eye(4);
			this.E_iw = eye(4);
			this.E_ip = eye(4);
			this.A0_ij = eye(6);
			this.phi = zeros(6,1);
			this.joint = {};
			countB = countB + 1;
			countCM = countCM + 1;
		end
		
		%%
		function setBodyTransform(this,E)
			% Sets the transform of this body wrt parent joint
			this.E0_ji = E;
			this.E0_ij = se3.inv(this.E0_ji);
			this.A0_ij = se3.Ad(this.E0_ij);
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
			this.E_wi = this.joint.E_wj*this.E0_ji;
			this.E_iw = se3.inv(this.E_wi);
			this.E_ip = eye(4);
			if ~isempty(this.joint.parent)
				this.E_ip = this.E_iw*this.joint.parent.body.E_wi;
			end
			% Body velocity
			this.phi = this.A0_ij*this.joint.V;
		end
		
		%%
		function [Mm,fm,Km,Dm] = computeMassGrav(this,grav,Mm,fm,Km,Dm)
			% Computes maximal mass matrix and force vector
			nm = redmax.Scene.countM();
			if nargout == 2
				if nargin == 2
					Mm = zeros(nm,nm);
					fm = zeros(nm,1);
				end
			else
				if nargin == 2
					Mm = zeros(nm,nm);
					fm = zeros(nm,1);
					Km = zeros(nm,nm);
					Dm = zeros(nm,nm);
				end
			end
			% Coriolis and gravity
			rows = this.idxM;
			M_i = diag(this.I_i);
			Mm(rows,rows) = M_i;
			adt = se3.ad(this.phi)';
			fcor = adt*M_i*this.phi;
			R_wi = this.E_wi(1:3,1:3);
			R_iw = R_wi';
			fgrav = zeros(6,1);
			mass = M_i(4,4);
			grav_i = R_iw*grav;
			fgrav(4:6) = mass*grav_i; % wrench in body space
			fm(rows) = fm(rows) + fcor + fgrav;
			if nargout == 2
				% Go to the next body
				if ~isempty(this.next)
					[Mm,fm] = this.next.computeMassGrav(grav,Mm,fm);
				end
			else
				% Derivative of gravity wrt position
				Km(rows(4:6),rows(1:3)) = Km(rows(4:6),rows(1:3)) + se3.brac(fgrav(4:6));
				% Derivative of Coriolis wrt velocity
				e1 = se3.brac([1 0 0]);
				e2 = se3.brac([0 1 0]);
				e3 = se3.brac([0 0 1]);
				z3 = [0 0 0]';
				Iw = this.I_i(1:3).*this.phi(1:3);
				mv = mass*this.phi(4:6);
				Dm(rows,rows) =  Dm(rows,rows) + adt*M_i - [
					e1*Iw, e2*Iw, e3*Iw, e1*mv, e2*mv, e3*mv
					e1*mv, e2*mv, e3*mv,    z3,    z3,    z3];
				% Go to the next body
				if ~isempty(this.next)
					[Mm,fm,Km,Dm] = this.next.computeMassGrav(grav,Mm,fm,Km,Dm);
				end
			end
		end
		
		%%
		function [fm,Km,Dm] = computeForce(this,fm,Km,Dm)
			% Computes maximal force vector and matrix
			nm = redmax.Scene.countM();
			if nargout == 1
				% Just the force
				if nargin == 1
					fm = zeros(nm,1);
				end
				% TODO: viscous damping
				% Go to the next body
				if ~isempty(this.next)
					fm = this.next.computeForce(fm);
				end
			else
				% Both force and derivatives
				if nargin == 1
					fm = zeros(nm,1);
					Km = zeros(nm);
					Dm = zeros(nm);
				end
				% TODO: viscous damping
				% Go to the next body
				if ~isempty(this.next)
					[fm,Km,Dm] = this.next.computeForce(fm,Km,Dm);
				end
			end
		end
		
		%%
		function [T,V] = computeEnergies(this,grav,T,V)
			% Computes kinetic and potential energies
			% The kinetic energy can also be computed by the joint as
			% T = 0.5*this.V'*this.I_j*this.V;
			T = T + 0.5*this.phi'*diag(this.I_i)*this.phi;
			V = V - this.I_i(6)*grav'*this.E_wi(1:3,4);
		end
		
		%%
		function draw(this)
			%s = this.getAxisSize();
			%if s > 0
			%	se3.drawAxis(this.E_wi,s);
			%end
			%text(this.E_wi(1,4),this.E_wi(2,4),this.E_wi(3,4),this.name);
			this.draw_();
			% Go to the next body
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
