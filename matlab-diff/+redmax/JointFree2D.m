classdef JointFree2D < redmax.Joint
	% 2D free joint in XY
	%
	
	%%
	properties
	end
	
	%%
	methods
		%%
		function this = JointFree2D(parent,body)
			this = this@redmax.Joint(parent,body,3);
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function update_(this,deriv)
			n = this.ndof;
			p = [this.q(1:2); 0];
			r = this.q(3);
			pdot = this.qdot(1:2);
			rdot = this.qdot(3);
			R = eye(3);
			s = sin(r);
			c = cos(r);
			R(1:2,1:2) = [c -s; s c];
			this.Q = [R p; 0 0 0 1];
			this.A = se3.Ad(this.Q);
			this.S(3,3) = 1;
			this.S(4:5,1:2) = [c s; -s c];
			
			pbrac = se3.brac(p);
			dRdq = zeros(3);
			dRdq(1:2,1:2) = [-s -c; c -s];
			
			pdotbrac = se3.brac([pdot;0]);
			Rdot = dRdq*rdot;
			this.Adot(1:3,1:3) = Rdot;
			this.Adot(4:6,4:6) = Rdot;
			this.Adot(4:6,1:3) = pdotbrac*R + pbrac*Rdot;
			this.Sdot(4:5,1:2) = [-s, c; -c, -s]*rdot;
			
			if deriv
				this.dSdotdq(4:5,1:2,3) = [-c, -s; s, -c]*rdot;
				
				e1brac = se3.brac([1 0 0]);
				e2brac = se3.brac([0 1 0]);
				this.dAdq(4:6,1:3,1) = e1brac*R;
				this.dAdq(4:6,1:3,2) = e2brac*R;
				this.dAdq(1:3,1:3,3) = dRdq;
				this.dAdq(4:6,1:3,3) = pbrac*dRdq;
				this.dAdq(4:6,4:6,3) = dRdq;
				
				this.dSdq = zeros(6,3,n);
				this.dSdq(4:5,1:2,3) = [-s c; -c -s];
				
				dRdotdq = zeros(3,3);
				dRdotdq(1:2,1:2) = [-c, s; -s, -c]*rdot;
				this.dAdotdq(4:6,1:3,1) = e1brac*Rdot;
				this.dAdotdq(4:6,1:3,2) = e2brac*Rdot;
				this.dAdotdq(1:3,1:3,3) = dRdotdq;
				this.dAdotdq(4:6,4:6,3) = dRdotdq;
				this.dAdotdq(4:6,1:3,3) = pdotbrac*dRdq + pbrac*dRdotdq;
			end
		end
	end
end
