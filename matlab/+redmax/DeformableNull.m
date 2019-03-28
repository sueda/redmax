classdef DeformableNull < redmax.Deformable
	% DeformableNull Does nothing
	
	%%
	methods
		%%
		function this = DeformableNull()
			this = this@redmax.Deformable();
			this.name = 'null';
		end
	end
end
