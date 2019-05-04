classdef Scene < handle
	%Scene Test scenes for redmax
	
	properties
		name
		bodies
		joints
		tspan
		hEuler
		drawHz
		grav
		waxis
	end
	
	methods
		function this = Scene()
			% Default values
			redmax.Scene.clear();
			this.name = '';
			this.bodies = {};
			this.joints = {};
			this.tspan = [0 2];
			this.hEuler = 1e-2;
			this.drawHz = 10;
			this.grav = [0 0 -980]';
			this.waxis = [];
		end
		
		function init(this)
			if usejava('jvm')
				colormap('default'); % Restore colormap mode
			end
			njoints = length(this.joints);
			joint0 = this.joints{1}; % Assumes joints{1} is root
			order = joint0.getTraversalOrder();
			this.joints = this.joints(order); % Reorders joints
			for i = njoints : -1 : 1
				% Leaf-to-root ordering is better for reducing matrix fill
				this.joints{i}.countDofs(); % This will also call Body.countDofs()
			end
			for i = 1 : njoints
				if i < njoints
					this.joints{i}.next = this.joints{i+1};
				end
				if i > 1
					this.joints{i}.prev = this.joints{i-1};
				end
			end
			joint0.update();
			nbodies = length(this.bodies);
			this.bodies = this.bodies(order); % same ordering as joints
			for i = 1 : nbodies
				this.bodies{i}.computeInertia();
				if i < nbodies
					this.bodies{i}.next = this.bodies{i+1}; %#ok<*SAGROW>
				end
			end
		end
		
		%%
		function draw(this,t)
			if t == 0
				clf;
				xlabel('X');
				ylabel('Y');
				zlabel('Z');
				axis equal;
				if ~isempty(this.waxis)
					axis(this.waxis);
				end
				ax = gca;
				ax.Clipping = 'off';
				grid on;
				view(3);
			end
			cla;
			hold on;
			this.bodies{1}.draw();
			this.joints{1}.draw();
			title(sprintf('t = %.4f', t));
			drawnow;
		end
	end
	
	%%
	methods (Static)
		%%
		function clear()
			global nm nr countCM CM;
			nm = 0;
			nr = 0;
			countCM = 1;
			if ~usejava('jvm')
				CM = zeros(1,3);
			else
				CM = colormap('lines');
			end
		end
		
		%%
		function out = countM(data)
			global nm;
			if nargin
				nm = data;
			end
			if isempty(nm)
				out = 0;
			else
				out = nm;
			end
		end
		
		%%
		function out = countR(data)
			global nr;
			if nargin
				nr = data;
			end
			if isempty(nr)
				out = 0;
			else
				out = nr;
			end
		end
	end
end
