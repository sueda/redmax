classdef JointTranslational < redmax.Joint
	% 3D translational joint
	
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
			n = 3;
			this.Q = eye(4);
			this.Q(1:3,4) = this.q;
			this.A = se3.Ad(this.Q);
			this.S = [zeros(3);eye(3)];
			
			this.Sdot = zeros(6,n);
			this.dSdotdq = zeros(6,n,n);
			
			this.dAdq = zeros(6,6,n);
			for k = 1 : 3
				ek = zeros(3,1);
				ek(k) = 1;
				ekbrac = se3.brac(ek);
				this.dAdq(4:6,1:3,k) = ekbrac;
			end
			this.dSdq = zeros(6,n,n);
			
			this.Adot = zeros(6,6);
			this.Adot(4:6,1:3) = se3.brac(this.qdot);
			this.dAdotdq = zeros(6,6,n);
		end
	end
end
