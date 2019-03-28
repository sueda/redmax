classdef JointUniversal < redmax.Joint
	% https://github.com/junggon/gear/blob/master/src/gjoint_universal.cpp
	% Assumes XY axes
	
	%%
	properties
	end
	
	%%
	methods
		%%
		function this = JointUniversal(parent,body)
			this = this@redmax.Joint(parent,body,2);
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function update_(this)
			q0 = this.q(1); q1 = this.q(2);
			dq0 = this.qdot(1); dq1 = this.qdot(2); %#ok<NASGU>
			c0 = cos(q0); c1 = cos(q1);
			s0 = sin(q0); s1 = sin(q1);
			this.Q = eye(4);
			this.Q(1:3,1) = [c1, s0*s1, -c0*s1]';
			this.Q(1:3,2) = [0, c0, s0]';
			this.Q(1:3,3) = [s1, -s0*c1, c0*c1]';
			%S[0] = c1; S[2] =  s1; S[7] = 1;
			this.S(1,1) = c1;
			this.S(3,1) = s1;
			this.S(2,2) = 1;
			%dS[0] = -s1*dq1; dS[2] = c1*dq1;
			this.Sdot(1,1) = -s1*dq1;
			this.Sdot(3,1) = c1*dq1;
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
