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
			this.Q = eye(4);
		end
	end
end
