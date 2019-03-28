classdef ConstraintNull < redmax.Constraint
	% ConstraintNull Does nothing
	
	%%
	methods
		%%
		function this = ConstraintNull()
			this = this@redmax.Constraint(0,0,0,0);
			this.name = 'null';
		end
	end
end
