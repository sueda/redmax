classdef JointTranslational < redmax.Joint
	% https://github.com/junggon/gear/blob/master/src/gjoint_translational.cpp
	
	%%
	properties
	end
	
	%%
	methods
		%%
		function this = JointTranslational(parent,body)
			this = this@redmax.Joint(parent,body,3);
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function update_(this)
			this.Q = eye(4);
			this.Q(1:3,4) = this.q(1:3);
			this.S(4,1) = 1.0;
			this.S(5,2) = 1.0;
			this.S(6,3) = 1.0;
		end
		
		%%
		function draw_(this)
			s = 2*min(this.body.sides);
			p = this.E_wj(1:3,4);
			x = this.E_wj(1:3,1);
			y = this.E_wj(1:3,2);
			z = this.E_wj(1:3,3);
			px0 = p - s*x;
			px1 = p + s*x;
			py0 = p - s*y;
			py1 = p + s*y;
			pz0 = p - s*z;
			pz1 = p + s*z;
			ps = [px0 px1 nan(3,1) py0 py1 nan(3,1) pz0 pz1];
			plot3(ps(1,:),ps(2,:),ps(3,:),'k','LineWidth',2);
		end
	end
end
