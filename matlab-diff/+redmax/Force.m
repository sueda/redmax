classdef Force < handle
	%Force Massless force acting on bodies
	
	%%
	properties
		next % Next force in traversal order
	end
	
	%%
	methods
		%%
		function this = Force()
		end
		
		%%
		function [fr,Kr,Dr,fm,Km,Dm] = computeValues(this,fr,Kr,Dr,fm,Km,Dm)
			nr = redmax.Scene.countR();
			nm = redmax.Scene.countM();
			if nargin == 1
				fr = zeros(nr,1);
				Kr = zeros(nr);
				Dr = zeros(nr);
				fm = zeros(nm,1);
				Km = zeros(nm);
				Dm = zeros(nm);
			end
			[fr,Kr,Dr,fm,Km,Dm] = this.computeValues_(fr,Kr,Dr,fm,Km,Dm);
			if ~isempty(this.next)
				[fr,Kr,Dr,fm,Km,Dm] = this.next.computeValues(fr,Kr,Dr,fm,Km,Dm);
			end
		end
		
		%%
		function V = computeEnergy(this,V)
			if nargin == 1
				V = 0;
			end
			V = this.computeEnergy_(V);
			if ~isempty(this.next)
				V = this.next.computeEnergy(V);
			end
		end
		
		%%
		function draw(this)
			this.draw_();
			if ~isempty(this.next)
				this.next.draw();
			end
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [fr,Kr,Dr,fm,Km,Dm] = computeValues_(this,fr,Kr,Dr,fm,Km,Dm) %#ok<*INUSL>
		end
		
		%%
		function V = computeEnergy_(this,V)
		end
		
		%%
		function draw_(this) %#ok<MANU>
		end
		
	end
end

