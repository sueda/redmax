classdef ForceNull < redmax.Force
	%ForceNull An empty force
	
	%%
	properties
	end
	
	methods
		function this = ForceNull()
			this = this@redmax.Force();
		end
	end
	
	%%
	methods (Access = protected)

	end
end
