classdef Scene < handle
	%Scene Test scenes for redmax
	
	properties
		name
		bodies
		joints
		forces
		tspan
		t % current time
		h % time step
		k % current step
		ts % times
		Vs % potential energies
		Ts % kinetic energies
		grav
		drawHz
		view
		waxis
		computeH
		plotH
		Hexpected
	end
	
	methods
		function this = Scene()
			% Default values
			redmax.Scene.clear();
			this.name = '';
			this.bodies = {};
			this.joints = {};
			this.forces = {};
			this.tspan = [0 1];
			this.h = 1e-2;
			this.t = 0;
			this.k = 0;
			this.ts = [];
			this.Vs = [];
			this.Ts = [];
			this.grav = [0 0 -980]';
			this.drawHz = 15;
			this.view = 3;
			this.waxis = [];
			this.computeH = true;
			this.plotH = true;
			this.Hexpected = zeros(1,2);
		end
		
		function init(this)
			if usejava('jvm')
				colormap('default'); % Restore colormap mode
			end
			njoints = length(this.joints);
			joint0 = this.joints{1}; % Assumes joints{1} is root

			% Leaf-to-root ordering is better for reducing matrix fill
			order = joint0.getTraversalOrder();
			this.joints = this.joints(order); % Reorders joints
			%for i = 1 : njoints % DEBUG ordering
			for i = njoints : -1 : 1
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
			
			nforces = length(this.forces);
			for i = 1 : nforces
				if i < nforces
					this.forces{i}.next = this.forces{i+1};
				end
			end
			if nforces == 0
				this.forces{1} = redmax.ForceNull();
			end
			
			% Set q0, qdot0, q1, qdot1
			for i = 1 : njoints
				this.joints{i}.q0 = this.joints{i}.q;
				this.joints{i}.qdot0 = this.joints{i}.qdot;
				this.joints{i}.q1 = this.joints{i}.q;
				this.joints{i}.qdot1 = this.joints{i}.qdot;
			end
			
			% Update all
			joint0.update();
			
			% Initialize bodies
			nbodies = length(this.bodies);
			this.bodies = this.bodies(order); % same ordering as joints
			for i = 1 : nbodies
				this.bodies{i}.computeInertia();
				if i < nbodies
					this.bodies{i}.next = this.bodies{i+1}; %#ok<*SAGROW>
				end
			end
			
			% Initial energy
			this.t = 0;
			this.k = 0;
			this.computeEnergies();
		end
		
		%%
		function computeEnergies(this)
			if this.computeH
				[T,V] = this.joints{1}.computeEnergies(this.grav);
				V = this.forces{1}.computeEnergy(V);
				this.Ts(end+1) = T;
				this.Vs(end+1) = V;
				this.ts(end+1) = this.t;
			end
		end
		
		%%
		function plotEnergies(this,itype)
			if this.computeH
				T = this.Ts;
				V = this.Vs - this.Vs(1);
				H = T + V;
				dH = H(end) - this.Hexpected(itype);
				if abs(dH) > 1e-2
					fprintf(2,'### FAIL: %.16e ###\n',H(end));
				else
					fprintf(1,'### PASS ###\n');
				end
				if this.plotH
					clf;
					plot(this.ts,T,'.-',this.ts,V,'.-',this.ts,H,'.-');
					a = axis();
					a(2) = this.ts(end);
					axis(a);
					grid on;
					legend('T','V','H');
				end
			end
		end
		
		%%
		function draw(this)
			if this.t == 0
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
				view(this.view); %#ok<CPROP>
			end
			timeToDraw = mod(this.t+1e-9,1/this.drawHz) < mod(this.t-this.h+2e-9,1/this.drawHz);
			if this.drawHz > 0 && timeToDraw
				cla;
				hold on;
				this.bodies{1}.draw();
				this.joints{1}.draw();
				title(sprintf('t = %.4f',this.t));
				drawnow;
			end
		end
		
		%%
		function test(this)
			% Test joint derivatives
			for i = 1 : length(this.joints)
				this.joints{i}.test();
			end
			
			fprintf('=== System ===\n');
			nm = redmax.Scene.countM();
			nr = redmax.Scene.countR();
			sqrteps = sqrt(eps);
			joint0 = this.joints{1};
			body0 = this.bodies{1};
			q = joint0.getQ();
			qdot = joint0.getQdot();
			[J,dJdq,Jdot,dJdotdq] = joint0.computeJacobian();
			[Mm,fm,Km,Dm] = body0.computeMassGrav(this.grav);
			[fm,Km,Dm] = body0.computeForce(fm,Km,Dm);
			[fr,Kr,Dr] = joint0.computeForce();
			
			% Inertia
			Mr = J'*Mm*J;
			dMrdq = zeros(nr,nr,nr);
			for i = 1 : nr
				tmp = J'*Mm*dJdq(:,:,i);
				dMrdq(:,:,i) = tmp' + tmp;
			end
			
			% Quadratic velocity vector
			fqvv = -J'*Mm*Jdot*qdot;
			Kqvv = zeros(nr,nr);
			Dqvv = -J'*Mm*Jdot;
			for i = 1 : nr
				dJdqk = dJdq(:,:,i);
				dJdotdqk = dJdotdq(:,:,i);
				Kqvv(:,i) = -dJdqk'*Mm*Jdot*qdot - J'*Mm*dJdotdqk*qdot;
				Dqvv(:,i) = Dqvv(:,i) - J'*Mm*dJdqk*qdot;
			end
			
			% Forces
			f = fr + J'*fm + fqvv;
			K = Kr + J'*Km*J + Kqvv;
			D = Dr + J'*Dm*J + Dqvv;
			for i = 1 : nr
				dJdqk = dJdq(:,:,i);
				K(:,i) = K(:,i) + dJdqk'*fm + J'*Dm*dJdqk*qdot;
			end
			
			%%%%%%%%%
			% TESTS %
			%%%%%%%%%
			
			% Test Jdot
			q_ = q + sqrteps*qdot;
			joint0.setQ(q_);
			joint0.update();
			J_ = joint0.computeJacobian();
			joint0.setQ(q);
			joint0.update();
			Jdot_ = (J_ - J)/sqrteps;
			redmax.Scene.printError('Jdot',Jdot_,Jdot);
			
			% Test dJdq and dJdotdq
			dJdq_ = zeros(nm,nr,nr);
			dJdotdq_ = zeros(nm,nr,nr);
			for i = 1 : nr
				q_ = q;
				q_(i) = q_(i) + sqrteps;
				joint0.setQ(q_);
				joint0.update();
				[J_,~,Jdot_] = joint0.computeJacobian();
				joint0.setQ(q);
				joint0.update();
				dJdq_(:,:,i) = (J_ - J)/sqrteps;
				dJdotdq_(:,:,i) = (Jdot_ - Jdot)/sqrteps;
			end
			redmax.Scene.printError('dJ/dq',dJdq_,dJdq);
			redmax.Scene.printError('dJdot/dq',dJdotdq_,dJdotdq);
			
			% Test inertia
			dMrdq_ = zeros(nr,nr,nr);
			for i = 1 : nr
				q_ = q;
				q_(i) = q_(i) + sqrteps;
				joint0.setQ(q_);
				joint0.update();
				J_ = joint0.computeJacobian();
				joint0.setQ(q);
				joint0.update();
				Mr_ = J_'*Mm*J_;
				dMrdq_(:,:,i) = (Mr_ - Mr)/sqrteps;
			end
			redmax.Scene.printError('dMr/dq',dMrdq_,dMrdq);
			
			% Test quadratic velocity vector
			Kqvv_ = zeros(nr,nr);
			Dqvv_ = zeros(nr,nr);
			for i = 1 : nr
				% Modify q
				q_ = q;
				q_(i) = q_(i) + sqrteps;
				joint0.setQ(q_);
				joint0.update();
				[J_,~,Jdot_] = joint0.computeJacobian();
				joint0.setQ(q);
				joint0.update();
				fqvv_ = -J_'*Mm*Jdot_*qdot;
				Kqvv_(:,i) = (fqvv_ - fqvv)/sqrteps;
				% Modify qdot
				qdot_ = qdot;
				qdot_(i) = qdot_(i) + sqrteps;
				joint0.setQdot(qdot_);
				joint0.update();
				[J_,~,Jdot_] = joint0.computeJacobian();
				joint0.setQdot(qdot);
				joint0.update();
				fqvv_ = -J_'*Mm*Jdot_*qdot_;
				Dqvv_(:,i) = (fqvv_ - fqvv)/sqrteps;
			end
			redmax.Scene.printError('Kqvv',Kqvv_,Kqvv);
			redmax.Scene.printError('Dqvv',Dqvv_,Dqvv);
			
			% Check K and D
			K_ = zeros(nr,nr);
			D_ = zeros(nr,nr);
			for i = 1 : nr
				% Modify q
				q_ = q;
				q_(i) = q_(i) + sqrteps;
				joint0.setQ(q_);
				joint0.update();
				[J_,~,Jdot_] = joint0.computeJacobian();
				[~,fm_] = body0.computeMassGrav(this.grav);
				fm_ = body0.computeForce(fm_);
				fr_ = joint0.computeForce();
				joint0.setQ(q);
				joint0.update();
				f_ = fr_ + J_'*fm_ - J_'*Mm*Jdot_*qdot;
				K_(:,i) = (f_ - f)/sqrteps;
				% Modify qdot
				qdot_ = qdot;
				qdot_(i) = qdot_(i) + sqrteps;
				joint0.setQdot(qdot_);
				joint0.update();
				[J_,~,Jdot_] = joint0.computeJacobian();
				[~,fm_] = body0.computeMassGrav(this.grav);
				fm_ = body0.computeForce(fm_);
				fr_ = joint0.computeForce();
				joint0.setQdot(qdot);
				joint0.update();
				f_ = fr_ + J_'*fm_ - J_'*Mm*Jdot_*qdot_;
				D_(:,i) = (f_ - f)/sqrteps;
			end
			redmax.Scene.printError('K',K_,K);
			redmax.Scene.printError('D',D_,D);
		end
	end
	
	%%
	methods (Static)
		%%
		function clear()
			global nm nr countCM CM countB;
			nm = 0;
			nr = 0;
			countCM = 1;
			countB = 1;
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
		function printError(name,v0,v1)
			thresh = 1e-6;
			s = size(v0);
			if length(s) < 3
				s = [s,1];
			end
			for k = 1 : s(3)
				e = norm(v1(:,:,k) - v0(:,:,k));
				nv0 = norm(v0(:,:,k));
				nv1 = norm(v1(:,:,k));
				if nv0 > 1e-4 && nv1 > 1e-4
					% Relative error
					e = e/min(nv0,nv1);
				end
				if s(3) == 1
					namek = name;
				else
					namek = sprintf('%s%d',name,k);
				end
				str = pad(namek,16);
				if abs(e) < thresh
					fprintf('%s%e\n',str,e);
				else
					fprintf(2,'%s%e\n',str,e);
				end
			end
		end
	end
end
