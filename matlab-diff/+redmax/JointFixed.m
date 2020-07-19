classdef JointFixed < redmax.Joint

	%%
	methods
		%%
		function this = JointFixed(parent,body)
			this = this@redmax.Joint(parent,body,0);
		end
	end
end
