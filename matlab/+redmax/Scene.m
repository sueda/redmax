classdef Scene < handle
	%Scene Test scenes for redmax
	
	properties
		name
		bodies
		joints
		constraints
		deformables
		springs
		tspan
		hEuler
		drawHz
		grav
		baumgarte
		fric
		mu % internal and external
		SPiterMax
		SPconv
		SPathresh
		SPreg
		waxis
		Hexpected
		sceneFcn
	end
	
	methods
		function this = Scene()
			% Default values
			redmax.Scene.clear();
			this.name = '';
			this.bodies = {};
			this.joints = {};
			this.constraints = {};
			this.deformables = {};
			this.springs = {};
			this.tspan = [0 2];
			this.hEuler = 1e-2;
			this.drawHz = 10;
			this.grav = [0 0 -980]';
			this.baumgarte = [5 5 5];
			this.fric = false;
			this.mu = [0.6 0.6];
			this.SPiterMax = 100;
			this.SPconv = 1e-3;
			this.SPathresh = 1e-10;
			this.SPreg = 1e-6;
			this.waxis = [];
			this.Hexpected = zeros(3,1);
			this.sceneFcn = @redmax.Scene.sceneFcnNull;
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
			joint0.generateContacts();
			nbodies = length(this.bodies);
			this.bodies = this.bodies(order); % same ordering as joints
			for i = 1 : nbodies
				this.bodies{i}.computeInertia();
				if i < nbodies
					this.bodies{i}.next = this.bodies{i+1}; %#ok<*SAGROW>
				end
			end
			ndeformables = length(this.deformables);
			for i = 1 : ndeformables
				this.deformables{i}.countDofs();
				this.deformables{i}.initGeometry();
				% Create attachment constraints
				this.constraints{end+1} = redmax.ConstraintAttachSpring(this.deformables{i});
				if i < ndeformables
					this.deformables{i}.next = this.deformables{i+1};
				end
			end
			if ndeformables == 0
				this.deformables{1} = redmax.DeformableNull();
			end
			nsprings = length(this.springs);
			for i = 1 : nsprings
				if i < nsprings
					this.springs{i}.next = this.springs{i+1};
				end
			end
			if nsprings == 0
				this.springs{1} = redmax.SpringNull();
			end
			nconstraints = length(this.constraints);
			for i = 1 : nconstraints
				this.constraints{i}.countDofs();
				if i < nconstraints
					this.constraints{i}.next = this.constraints{i+1};
				end
			end
			if nconstraints == 0
				this.constraints{1} = redmax.ConstraintNull();
			end
			% Create static contacts for (loop-closing) joints
			this.constraints{1}.generateContactsJoint();
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
			this.sceneFcn(t,this);
			cla;
			hold on;
			this.bodies{1}.draw();
			this.joints{1}.draw();
			this.deformables{1}.draw();
			this.springs{1}.draw();
			this.constraints{1}.draw();
			title(sprintf('t = %.4f', t));
			drawnow;
		end
	end
	
	methods (Static)
		function sceneFcnNull(scene,t) %#ok<INUSD>
		end
	end
	
	%%
	methods (Static)
		%%
		function clear()
			global nm nr nt nem ner nim nir countB countD countCM CM;
			nm = 0;
			nr = 0;
			nt = 0;
			nem = 0;
			ner = 0;
			nim = 0;
			nir = 0;
			countB = 1;
			countD = 1;
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
		
		%%
		function out = countT(data)
			global nt;
			if nargin
				nt = data;
			end
			if isempty(nt)
				out = 0;
			else
				out = nt;
			end
		end
		
		%%
		function out = countEM(data)
			global nem;
			if nargin
				nem = data;
			end
			if isempty(nem)
				out = 0;
			else
				out = nem;
			end
		end
		
		%%
		function out = countER(data)
			global ner;
			if nargin
				ner = data;
			end
			if isempty(ner)
				out = 0;
			else
				out = ner;
			end
		end
		
		%%
		function out = countIM(data)
			global nim;
			if nargin
				nim = data;
			end
			if isempty(nim)
				out = 0;
			else
				out = nim;
			end
		end
		
		%%
		function out = countIR(data)
			global nir;
			if nargin
				nir = data;
			end
			if isempty(nir)
				out = 0;
			else
				out = nir;
			end
		end
	end
end
