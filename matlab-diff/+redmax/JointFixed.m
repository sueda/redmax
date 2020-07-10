classdef JointFixed < redmax.Joint

	%%
	methods
		%%
		function this = JointFixed(parent,body)
			this = this@redmax.Joint(parent,body,0);
		end
	end
		
	%%
	methods (Access = protected)
		%%
		function update_(this)
			this.Q = eye(4,4);
			this.A = eye(6,6);
			this.dAdq = zeros(6,6,0);
			this.Adot = zeros(6,6);
			this.dAdotdq = zeros(6,6,0);
			this.S = zeros(6,0);
			this.dSdq = zeros(6,0,0);
			this.Sdot = zeros(6,0);
			this.dSdotdq = zeros(6,0,0);
		end
	end
end
