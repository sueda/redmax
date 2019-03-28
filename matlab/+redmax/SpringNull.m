classdef SpringNull < redmax.Spring
	%SpringNull Empty spring, created if there are no other springs

	%%
	methods
		function this = SpringNull()
			this = this@redmax.Spring();
		end
	end
end
