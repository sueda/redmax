classdef (Abstract) Force < handle
	%Force Massless force acting on bodies
	
	%%
	properties
		name % optional name
		next % Next force in traversal order
	end
	
	%%
	methods
		%%
		function this = Force()
		end
		
		%%
		function init(this)
			this.init_();
			% Go to the next force
			if ~isempty(this.next)
				this.next.init();
			end
		end
		
		%%
		function [fr,fm,Kr,Km,Dr,Dm] = computeValues(this,fr,fm,Kr,Km,Dr,Dm)
			nr = redmax.Scene.countR();
			nm = redmax.Scene.countM();
			if nargout == 2
				% Just the forces
				if nargin == 1
					fr = zeros(nr,1);
					fm = zeros(nm,1);
				end
				[fr,fm] = this.computeValues_(fr,fm);
				% Go to the next force
				if ~isempty(this.next)
					[fr,fm] = this.next.computeValues(fr,fm);
				end
			else
				% Forces and derivatives
				if nargin == 1
					fr = zeros(nr,1);
					fm = zeros(nm,1);
					Kr = zeros(nr);
					Km = zeros(nm);
					Dr = zeros(nr);
					Dm = zeros(nm);
				end
				[fr,fm,Kr,Km,Dr,Dm] = this.computeValues_(fr,fm,Kr,Km,Dr,Dm);
				% Go to the next force
				if ~isempty(this.next)
					[fr,fm,Kr,Km,Dr,Dm] = this.next.computeValues(fr,fm,Kr,Km,Dr,Dm);
				end
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
		function init_(this) %#ok<MANU>
		end
		
		%%
		function [fr,fm,Kr,Km,Dr,Dm] = computeValues_(this,fr,fm,Kr,Km,Dr,Dm) %#ok<*INUSL>
		end
		
		%%
		function V = computeEnergy_(this,V)
		end
		
		%%
		function draw_(this) %#ok<MANU>
		end
		
	end
end

