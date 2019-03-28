classdef JointPlanar < redmax.Joint
	% https://github.com/junggon/gear/blob/master/src/gjoint_planar.cpp
	% Assumes X-Y plane
	
	%%
	properties
	end
	
	%%
	methods
		%%
		function this = JointPlanar(parent,body)
			this = this@redmax.Joint(parent,body,2);
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function update_(this)
			this.Q = eye(4);
			this.Q(1:2,4) = this.q(1:2);
			this.S(4,1) = 1.0;
			this.S(5,2) = 1.0;
		end
		
		%%
		function draw_(this)
			s = 2*min(this.body.sides);
			p = this.E_wj(1:3,4);
			x = this.E_wj(1:3,1);
			y = this.E_wj(1:3,2);
			px0 = p - s*x;
			px1 = p + s*x;
			py0 = p - s*y;
			py1 = p + s*y;
			ps = [px0 px1 nan(3,1) py0 py1];
			plot3(ps(1,:),ps(2,:),ps(3,:),'k','LineWidth',2);
		end
	end
end
