classdef JointFree3D < redmax.Joint
	% 3D free joint
	% This composite joint uses JointTranslational as S1 and
	% JointSphericalExp as S2, unlike JointFree. This version behaves
	% better with Euler.
	%
	
	%%
	properties
		joint1
		joint2
	end
	
	%%
	methods
		%%
		function this = JointFree3D(parent,body)
			this = this@redmax.Joint(parent,body,6);
			this.joint1 = redmax.JointTranslational([],body);
			this.joint2 = redmax.JointSphericalExp([],body);
			body.joint = this;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function reparam_(this)
			this.joint2.reparam_();
			this.q(4:6) = this.joint2.q;
			this.qdot(4:6) = this.joint2.qdot;
		end
		
		%%
		function update_(this)
			this.joint1.q = this.q(1:3);
			this.joint2.q = this.q(4:6);
			this.joint1.qdot = this.qdot(1:3);
			this.joint2.qdot = this.qdot(4:6);
			this.joint1.update_();
			this.joint2.update_();
			
			Q1 = this.joint1.Q;
			Q2 = this.joint2.Q;
			this.Q = Q1*Q2;
			
			S2 = this.joint2.S;
			S2hat = S2(1:3,1:3);
			R_21 = Q2(1:3,1:3)';
			this.S(4:6,1:3) = R_21;
			this.S(1:3,4:6) = S2hat;
			
			dS2 = this.joint2.Sdot;
			this.Sdot(4:6,1:3) = -se3.brac(S2hat*this.qdot(4:6))*R_21;
			this.Sdot(1:3,4:6) = dS2(1:3,1:3);
		end
	end
end
