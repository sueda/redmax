classdef (Abstract) Deformable < handle
	
	%%
	properties
		name
		color
		damping
		next
	end
	
	%%
	methods
		%%
		function this = Deformable()
			global countD countCM CM;
			this.name = ['deformable',num2str(countD)];
			this.color = CM(mod(countCM-1,size(CM,1))+1,:);
			this.damping = 0;
			countD = countD + 1;
			countCM = countCM + 1;
		end
		
		%%
		function setDamping(this,damping)
			% Sets viscous damping
			this.damping = damping;
		end
		
		%%
		function initGeometry(this)
			this.initGeometry_();
			if ~isempty(this.next)
				this.next.initGeometry();
			end
		end
		
		%%
		function countDofs(this)
			this.countDofs_();
			if ~isempty(this.next)
				this.next.countDofs();
			end
		end
		
		%%
		function y = gatherDofs(this,y)
			y = this.gatherDofs_(y);
			if ~isempty(this.next)
				y = this.next.gatherDofs(y);
			end
		end
		
		%%
		function ydot = gatherDDofs(this,ydot)
			ydot = this.gatherDDofs(ydot);
			if ~isempty(this.next)
				ydot = this.next.gatherDDofs(ydot);
			end
		end
		
		%%
		function scatterDofs(this,y)
			this.scatterDofs_(y);
			if ~isempty(this.next)
				this.next.scatterDofs(y);
			end
		end
		
		%%
		function scatterDDofs(this,ydot)
			this.scatterDDofs_(ydot);
			if ~isempty(this.next)
				this.next.scatterDDofs(ydot);
			end
		end
		
		%%
		function [J,Jdot] = computeJacobian(this,J,Jdot)
			[J,Jdot] = this.computeJacobian_(J,Jdot);
			if ~isempty(this.next)
				[J,Jdot] = this.next.computeJacobian(J,Jdot);
			end
		end
		
		%%
		function [M,f] = computeMassGrav(this,grav,M,f)
			[M,f] = this.computeMassGrav_(grav,M,f);
			if ~isempty(this.next)
				[M,f] = this.next.computeMassGrav(grav,M,f);
			end
		end
		
		%%
		function [f,D] = computeForceDamping(this,f,D)
			[f,D] = this.computeForceDamping_(f,D);
			if ~isempty(this.next)
				[f,D] = this.next.computeForceDamping(f,D);
			end
		end
		
		%%
		function [T,V] = computeEnergies(this,grav,T,V)
			[T,V] = this.computeEnergies_(grav,T,V);
			if ~isempty(this.next)
				[T,V] = this.next.computeEnergies(grav,T,V);
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
		function initGeometry_(this) %#ok<MANU>
		end
		
		%%
		function countDofs_(this) %#ok<MANU>
		end
		
		%%
		function y = gatherDofs_(this,y) %#ok<INUSL>
		end
		
		%%
		function ydot = gatherDDofs_(this,ydot) %#ok<INUSL>
		end
		
		%%
		function scatterDofs_(this,y) %#ok<INUSD>
		end
		
		%%
		function scatterDDofs_(this,ydot) %#ok<INUSD>
		end
		
		%%
		function [J,Jdot] = computeJacobian_(this,J,Jdot) %#ok<INUSL>
		end
		
		%%
		function [M,f] = computeMassGrav_(this,grav,M,f) %#ok<INUSL>
		end
		
		%%
		function [f,D] = computeForceDamping_(this,f,D) %#ok<INUSL>
		end
		
		%%
		function [T,V] = computeEnergies_(this,grav,T,V) %#ok<INUSL>
		end
		
		%%
		function draw_(this) %#ok<MANU>
		end
	end
end

