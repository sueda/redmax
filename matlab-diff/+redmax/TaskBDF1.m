classdef (Abstract) TaskBDF1 < handle
	%TaskBDF1 A user-defined objective to be minimized
	
	%%
	properties
		scene
		p % parameters
		P % objective value
		wreg % regularizer weight
		dgdp % 
		dPdq %
	end
	
	%%
	methods
		%%
		function this = TaskBDF1(scene,nparams)
			this.scene = scene;
			this.p = zeros(nparams,1);
			this.P = 0;
			this.wreg = 1;
			this.dgdp = [];
			this.dPdq = {};
		end
		
		%% Initializes the task
		function init(this)
			nr = redmax.Scene.countR();
			np = length(this.p);
			this.P = 0;
			nsteps = this.scene.nsteps;
			this.dPdq = cell(1,nsteps);
			this.dgdp = sparse(nsteps*nr,np);
		end
		
		%% Applies the parameters to the current step
		function applyStep(this) %#ok<MANU>
		end
		
		%% Computes objective value and derivatives after a sim step
		function calcStep(this) %#ok<MANU>
		end
		
		%% Computes objective value and derivatives after the sim
		function [P,dPdp] = calcFinal(this)
			nr = redmax.Scene.countR();
			
			% Assume that all tasks have a simple regularizer
			P = this.P + this.wreg*0.5*(this.p'*this.p);
			
			% Compute derivative
			nsteps = this.scene.nsteps;
			z = zeros(nsteps*nr,1);
			h = this.scene.h;
			for k = this.scene.nsteps : -1 : 1
				yk = this.dPdq{k};
				kk0 = (k-1)*nr + (1:nr); % indices for time step k
				kk1 = kk0 + nr;          % indices for time step k+1
				kk2 = kk1 + nr;          % indices for time step k+2
				if k < nsteps
					% Add contributions from step k+1
					M = this.scene.history(k+1).M;
					D = this.scene.history(k+1).D;
					block = -2*M + h*D;
					yk = yk - block'*z(kk1);
				end
				if k < nsteps - 1
					% Add contributions from step k+2
					M = this.scene.history(k+2).M;
					block = M;
					yk = yk - block'*z(kk2);
				end
				% Solve by the diagonal
				Hl = this.scene.history(k).Hl;
				Hu = this.scene.history(k).Hu;
				Hp = this.scene.history(k).Hp;
				zkk0(Hp) = Hl'\(Hu'\yk);
				z(kk0) = zkk0;
			end
 			dPdp = this.wreg*this.p' - z'*this.dgdp;
		end
		
		%% Draws the current task state
		function draw(this) %#ok<MANU>
		end
	end
end
