classdef ForceCable < redmax.ForceSpringMultiPointGeneric
	%ForceCable A cable routed through multiple points
	
	%%
	properties
		L
		stiffness
		damping
	end
	
	methods
		function this = ForceCable()
			this = this@redmax.ForceSpringMultiPointGeneric();
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
			this.L = 0;
			npts = length(this.bodies);
			for k = 1 : npts - 1
				E1 = eye(4);
				E2 = eye(4);
				xl1 = this.xls{k};
				xl2 = this.xls{k+1};
				if ~isempty(this.bodies{k})
					E1 = this.bodies{k}.E_wi;
				end
				if ~isempty(this.bodies{k+1})
					E2 = this.bodies{k+1}.E_wi;
				end
				R1 = E1(1:3,1:3);
				R2 = E2(1:3,1:3);
				p1 = E1(1:3,4);
				p2 = E2(1:3,4);
				xw1 = R1*xl1 + p1;
				xw2 = R2*xl2 + p2;
				dx = xw2 - xw1;
				this.L = this.L + norm(dx);
			end
		end
		
		%%
		function [V,f,dfdl,dfdldot] = computeSpringForce(this,l,ldot)
			% If the force is positive, the spring will try to contract
			strain = (l - this.L)/this.L;
			dstrain = ldot/this.L;
			if strain > 0
				V = (this.stiffness/2)*strain^2*this.L;
				f = this.stiffness*strain + this.damping*dstrain;
				dfdl = this.stiffness/this.L;
				dfdldot = this.damping/this.L;
			else
				V = 0;
				f = 0;
				dfdl = 0;
				dfdldot = 0;
			end
		end
	end
end
