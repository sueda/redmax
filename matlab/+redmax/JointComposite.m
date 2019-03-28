classdef JointComposite < redmax.Joint
	% https://github.com/junggon/gear/blob/master/src/gjoint_composite.cpp
	
	%%
	properties
		joint1
		joint2
	end
	
	%%
	methods
		%%
		function this = JointComposite(parent,body,joint1,joint2)
			this = this@redmax.Joint(parent,body,joint1.ndof+joint2.ndof);
			this.joint1 = joint1;
			this.joint2 = joint2;
			body.joint = this;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function reparam_(this)
			ndof1 = this.joint1.ndof;
			ndof2 = this.joint2.ndof;
			idx1 = 1:ndof1;
			idx2 = ndof1+1:ndof1+ndof2;
			this.joint1.reparam_();
			this.joint2.reparam_();
			this.q(idx1) = this.joint1.q;
			this.q(idx2) = this.joint2.q;
			this.qdot(idx1) = this.joint1.qdot;
			this.qdot(idx2) = this.joint2.qdot;
		end
		
		%%
		function update_(this)
			ndof1 = this.joint1.ndof;
			ndof2 = this.joint2.ndof;
			idx1 = 1:ndof1;
			idx2 = ndof1+1:ndof1+ndof2;
			this.joint1.q = this.q(idx1);
			this.joint2.q = this.q(idx2);
			this.joint1.qdot = this.qdot(idx1);
			this.joint2.qdot = this.qdot(idx2);
			this.joint1.update_();
			this.joint2.update_();
			
			Q1 = this.joint1.Q;
			Q2 = this.joint2.Q;
			this.Q = Q1*Q2;
			
			S1 = this.joint1.S;
			S2 = this.joint2.S;
			Ad_21 = se3.Ad(Q2);
			this.S(1:6,idx1) = Ad_21*S1;
			this.S(1:6,idx2) = S2;
			
			dS1 = this.joint1.Sdot;
			dS2 = this.joint2.Sdot;
			Sdq2 = S2*this.joint2.qdot;
			this.Sdot(1:6,idx1) = -se3.ad(Sdq2)*Ad_21*S1 + Ad_21*dS1;
			this.Sdot(1:6,idx2) = dS2;
		end
		
		%%
		function draw_(this) %#ok<MANU>
			% TODO
		end
	end
end
