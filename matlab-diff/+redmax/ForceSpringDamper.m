classdef ForceSpringDamper < redmax.ForceSpringGeneric
	%ForceSpringDamper A damped spring force between two points
	
	%%
	properties
		L
		stiffness
		damping
	end
	
	methods
		function this = ForceSpringDamper(body1,x_1,body2,x_2)
			this = this@redmax.ForceSpringGeneric(body1,x_1,body2,x_2);
			this.stiffness = 1;
			this.damping = 1;
			this.L = 0;
		end
		
		%%
		function setStiffness(this,stiffness)
			this.stiffness = stiffness;
		end
		
		%%
		function setDamping(this,damping)
			this.damping = damping;
		end
		
		%%
		function setRetLength(this,L)
			this.L = L;
		end
	end
	
	methods (Access = protected)
		%%
		function init_(this)
			% Compute rest length, unless it has been set manually
			if this.L > 0
				return;
			end
			E1 = eye(4);
			E2 = eye(4);
			xl1 = this.x_1;
			xl2 = this.x_2;
			if ~isempty(this.body1)
				E1 = this.body1.E_wi;
			end
			if ~isempty(this.body2)
				E2 = this.body2.E_wi;
			end
			R1 = E1(1:3,1:3);
			R2 = E2(1:3,1:3);
			p1 = E1(1:3,4);
			p2 = E2(1:3,4);
			xw1 = R1*xl1 + p1;
			xw2 = R2*xl2 + p2;
			dx = xw2 - xw1;
			this.L = norm(dx);
		end
		
		%%
		function [V,f,dfdl,dfdldot] = computeSpringForce(this,l,ldot)
			% If the force is positive, the spring will try to contract
			strain = (l-this.L)/this.L;
			dstrain = ldot/this.L;
			V = (this.stiffness/2)*strain^2*this.L;
			f = this.stiffness*strain + this.damping*dstrain;
			dfdl = this.stiffness/this.L;
			dfdldot = this.damping/this.L;
		end
	end
end
