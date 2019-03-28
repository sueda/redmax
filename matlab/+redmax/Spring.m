classdef Spring < handle
	%Spring Mass-less spring acting on bodies
	
	%%
	properties
		next % Next spring in traversal order
	end
	
	%%
	methods
		%%
		function this = Spring()
		end
		
		%%
		function [f,K,D] = computeForceStiffnessDamping(this,f,K,D)
			nm = redmax.Scene.countM();
			if nargin == 1
				f = zeros(nm,1);
				K = zeros(nm);
				D = zeros(nm);
			elseif nargin == 2
				K = zeros(nm);
				D = zeros(nm);
			end
			if isempty(K)
				K = zeros(nm);
			end
			[f,K,D] = this.computeForceStiffnessDamping_(f,K,D);
			if ~isempty(this.next)
				[f,K,D] = this.next.computeForceStiffnessDamping(f,K,D);
			end
		end
		
		%%
		function V = computeEnergy(this,V)
			V = this.computeEnergy_(V);
			if ~isempty(this.next)
				V = this.next.computeEnergy(V);
			end
		end
		
		%%
		function y = computeStiffnessProd(this,x,y)
			% Computes y = K*x
			if nargin == 2
				nm = redmax.Scene.countM();
				y = zeros(nm,1);
			end
			y = this.computeStiffnessProd_(x,y);
			if ~isempty(this.next)
				y = this.next.computeStiffnessProd(x,y);
			end
		end
		
		%%
		function y = computeDampingProd(this,x,y)
			% Computes y = D*x
			if nargin == 2
				nm = redmax.Scene.countM();
				y = zeros(nm,1);
			end
			y = this.computeDampingProd_(x,y);
			if ~isempty(this.next)
				y = this.next.computeDampingProd(x,y);
			end
		end
		
		%%
		function u = getStiffnessRank1(this)
			nm = redmax.Scene.countM();
			u = zeros(nm,1);
			u = this.getStiffnessRank1_(u);
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
		function [f,K,D] = computeForceStiffnessDamping_(this,f,K,D) %#ok<*INUSL>
		end
		
		%%
		function V = computeEnergy_(this,V)
		end
		
		%%
		function y = computeStiffnessProd_(this,x,y)
		end
		
		%%
		function y = computeDampingProd_(this,x,y)
		end
		
		%%
		function u = getStiffnessRank1_(this,u)
		end
		
		%%
		function draw_(this) %#ok<MANU>
		end
		
	end
end

