classdef JointPlanar < redmax.Joint
	
	%%
	properties
		plane
	end
	
	%%
	methods
		%%
		function this = JointPlanar(parent,body,plane)
			this = this@redmax.Joint(parent,body,2);
			if nargin == 2
				plane = [1 0 0; 0 1 0]';
			end
			this.plane(:,1) = plane(:,1)/norm(plane(:,1));
			this.plane(:,2) = plane(:,2)/norm(plane(:,2));
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function update_(this)
			n = 2;
			B = this.plane;
			this.Q = eye(4);
			this.Q(1:3,4) = B*this.q;
			this.A = se3.Ad(this.Q);
			this.S = [zeros(3,2);B];
			
			this.Sdot = zeros(6,n);
			this.dSdotdq = zeros(6,n,n);
			
			this.dAdq = zeros(6,6,n);
			for k = 1 : 2
				bkbrac = se3.brac(B(:,k));
				this.dAdq(4:6,1:3,k) = bkbrac;
			end
			this.dSdq = zeros(6,n,n);
			
			this.Adot = zeros(6,6);
			this.Adot(4:6,1:3) = se3.brac(B*this.qdot);
			this.dAdotdq = zeros(6,6,n);
		end
	end
end
