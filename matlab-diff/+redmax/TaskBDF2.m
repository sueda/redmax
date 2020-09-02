classdef (Abstract) TaskBDF2 < handle
	%TaskBDF2 A user-defined objective to be minimized
	
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
		function this = TaskBDF2(scene,nparams)
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
			% We are not computing dg/dqa, assuming that dP/dqa is 0. In
			% other words, we're assuming that the objective does not
			% depend on time step alpha.
			nsteps = this.scene.nsteps;
			z = zeros(nsteps*nr,1);
			h = this.scene.h;
			a = (2-sqrt(2))/2;
			for k = this.scene.nsteps : -1 : 1
				yk = this.dPdq{k};
				kk0 = (k-1)*nr + (1:nr); % indices for time step k
				kk1 = kk0 + nr;          % indices for time step k+1
				kk2 = kk1 + nr;          % indices for time step k+2
				kk3 = kk2 + nr;          % indices for time step k+3
				kk4 = kk3 + nr;          % indices for time step k+4
				if k < nsteps
					% Add contributions from step k+1
					M = this.scene.history(k+1).M;
					D = this.scene.history(k+1).D;
					if k == 1
						block = -((8/(9*a)) + (4/3))*M + (8/9)*h*D;
					else
						block = -(8/3)*M + (8/9)*h*D;
					end
					yk = yk - block'*z(kk1);
				end
				if k < nsteps - 1
					% Add contributions from step k+2
					M = this.scene.history(k+2).M;
					D = this.scene.history(k+2).D;
					if k == 1
						block = ((2/(9*a)) + (19/9))*M - (2/9)*h*D;
					else
						block = (22/9)*M - (2/9)*h*D;
					end
					yk = yk - block'*z(kk2);
				end
				if k < nsteps - 2
					% Add contributions from step k+3
					M = this.scene.history(k+3).M;
					block = -(8/9)*M;
					yk = yk - block'*z(kk3);
				end
				if k < nsteps - 3
					% Add contributions from step k+4
					M = this.scene.history(k+4).M;
					block = (1/9)*M;
					yk = yk - block'*z(kk4);
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
