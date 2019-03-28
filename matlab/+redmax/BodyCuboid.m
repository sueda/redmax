classdef BodyCuboid < redmax.Body
	%%
	properties
		sides    % Side lengths
	end
	
	%%
	methods
		%%
		function this = BodyCuboid(density,sides)
			this = this@redmax.Body(density);
			this.sides = sides;
		end
		
		%%
		function computeInertia_(this)
			% Computes inertia at body
			this.I_i = se3.inertiaCuboid(this.sides,this.density);
		end
		
		%%
		function draw_(this)
			[F,V] = se3.patchCuboid(this.E_wi,this.sides);
			patch('Faces',F,'Vertices',V,'FaceColor',this.color);
		end
		
		%%
		function s = getAxisSize(this)
			s = min(this.sides);
		end
	end
end
