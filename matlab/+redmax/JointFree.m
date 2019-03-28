classdef JointFree < redmax.Joint
	% https://github.com/junggon/gear/blob/master/src/gjoint_free.cpp
	% Copied from GJointFreeST
	%
	
	%%
	properties
		jointS
		jointT
	end
	
	%%
	methods
		%%
		function this = JointFree(parent,body)
			this = this@redmax.Joint(parent,body,6);
			this.jointS = redmax.JointSphericalExp([],body);
			this.jointT = redmax.JointTranslational([],body);
			body.joint = this;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function reparam_(this)
			this.jointS.reparam_();
			this.q(1:3) = this.jointS.q;
			this.qdot(1:3) = this.jointS.qdot;
		end
		
		%%
		function update_(this)
			this.jointS.q = this.q(1:3);
			this.jointT.q = this.q(4:6);
			this.jointS.qdot = this.qdot(1:3);
			this.jointT.qdot = this.qdot(4:6);
			this.jointS.update_();
			this.jointT.update_();
			
			% T1 = spherical_joint.T = SE3(R, 0)
			% T2 = translational_joint.T = SE3(eye, p)
			% T = T1*T2 = SE3(R, R*p)
			Q1 = this.jointS.Q;
			Q2 = this.jointT.Q;
			this.Q = Q1*Q2;
			
			% S.Push(0, 0, S1w);
			% S.Push(3, 0, -Cross(p, S1w));
			% S.Push(0, 3, translational_joint.S);
			S1 = this.jointS.S;
			S2 = this.jointT.S;
			S1w = S1(1:3,1:3);
			p = Q2(1:3,4);
			pbrac = se3.brac(p);
			this.S(1:3,1:3) = S1w;
			this.S(4:6,1:3) = -pbrac*S1w;
			this.S(1:6,4:6) = S2;
			
			% dS.Push(0, 0, dS1w);
			% dS.Push(3, 0, -(Cross(S2dq2v, S1w) + Cross(p, dS1w)));
			% dS.Push(0, 3, translational_joint.dS);
			dS1 = this.jointS.Sdot;
			dS2 = this.jointT.Sdot;
			dS1w = dS1(1:3,1:3);
			this.Sdot(1:3,1:3) = dS1w;
			this.Sdot(4:6,1:3) = -(se3.brac(this.qdot(4:6))*S1w + pbrac*dS1w);
			this.Sdot(1:6,4:6) = dS2;
		end
	end
end
