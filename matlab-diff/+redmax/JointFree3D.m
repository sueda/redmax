classdef JointFree3D < redmax.Joint
	% 3D free joint
	% This composite joint uses JointTranslational as S1 and
	% JointSpherical as S2.
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
			this.joint2 = redmax.JointSpherical([],body);
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
			
			% Extract from translational joint
			p = this.joint1.q;
			pdot = this.joint1.qdot;
			
			% Extract from spherical joint
			rdot = this.joint2.qdot;
			R = this.joint2.Q(1:3,1:3);
			Rdot = this.joint2.Adot(1:3,1:3);
			T = this.joint2.S(1:3,1:3);
			Tdot = this.joint2.Sdot(1:3,1:3);
			dRdr = zeros(3,3,3);
			dRdotdr = zeros(3,3,3);
			dTdr = zeros(3,3,3);
			dTdotdr = zeros(3,3,3);
			for k = 1 : 3
				dRdr(:,:,k) = this.joint2.dAdq(1:3,1:3,k);
				dRdotdr(:,:,k) = this.joint2.dAdotdq(1:3,1:3,k);
				dTdr(:,:,k) = this.joint2.dSdq(1:3,:,k);
				dTdotdr(:,:,k) = this.joint2.dSdotdq(1:3,:,k);
			end
			
			n = 6;
			this.Q = eye(4,4);
			this.A = zeros(6,6);
			this.dAdq = zeros(6,6,n);
			this.Adot = zeros(6,6);
			this.dAdotdq = zeros(6,6,n);
			this.S = zeros(6,n);
			this.dSdq = zeros(6,n,n);
			this.Sdot = zeros(6,n);
			this.dSdotdq = zeros(6,n,n);
			
			% Q and A
			this.Q(1:3,1:3) = R;
			this.Q(1:3,4) = p;
			this.A(1:3,1:3) = R;
			this.A(4:6,4:6) = R;
			pbrac = se3.brac(p);
			this.A(4:6,1:3) = pbrac*R;
			
			% dAdq
			for k = 1 : 3
				ek = zeros(3,1);
				ek(k) = 1;
				ekbrac = se3.brac(ek);
				dRdrk = dRdr(:,:,k);
				this.dAdq(1:3,1:3,3+k) = dRdrk;
				this.dAdq(4:6,4:6,3+k) = dRdrk;
				this.dAdq(4:6,1:3,k) = ekbrac*R;
				this.dAdq(4:6,1:3,3+k) = pbrac*dRdrk;
			end
			
			% Adot
			pdotbrac = se3.brac(pdot);
			this.Adot(1:3,1:3) = Rdot;
			this.Adot(4:6,4:6) = Rdot;
			this.Adot(4:6,1:3) = pdotbrac*R + pbrac*Rdot;
			
			% dAdotdq
			for k = 1 : 3
				ek = zeros(3,1);
				ek(k) = 1;
				ekbrac = se3.brac(ek);
				dRdotdrk = dRdotdr(:,:,k);
				this.dAdotdq(4:6,1:3,k) = ekbrac*Rdot;
				this.dAdotdq(1:3,1:3,3+k) = dRdotdrk;
				this.dAdotdq(4:6,4:6,3+k) = dRdotdrk;
				this.dAdotdq(4:6,1:3,3+k) = pdotbrac*dRdr(:,:,k) + pbrac*dRdotdrk;
			end
			
			% S
			this.S(4:6,1:3) = R';
			this.S(1:3,4:6) = T;
			
			% dSdq
			for k = 1 : 3
				dRdrk = dRdr(:,:,k);
				this.dSdq(4:6,1:3,3+k) = dRdrk';
				this.dSdq(1:3,4:6,3+k) = dTdr(:,:,k);
			end
			
			% Sdot
			this.Sdot(4:6,1:3) = Rdot';
			this.Sdot(1:3,4:6) = Tdot;
			
			% dSdotdq
			tmp = se3.brac(T*rdot);
			for k = 1 : 3
				this.dSdotdq(4:6,1:3,3+k) = -se3.brac(dTdr(:,:,k)*rdot)*R' - tmp*dRdr(:,:,k)';
				this.dSdotdq(1:3,4:6,3+k) = dTdotdr(:,:,k);
			end
		end
	end
end
