classdef Scene < handle
	%Scene Test scenes for redmax
	
	properties
		name % scene name
		bodies % list of bodies
		joints % list of joints
		forces % list of forces
		tEnd % end time
		qInit % initial positions
		qdotInit % initial velocities
		t % current time
		h % time step
		k % current step
		T0 % initial kinetic energy
		V0 % initial potential energy
		history % step history: does not include the initial state
		nsteps % number of steps to take
		grav % gravity
		drawHz % refresh rate (0 for no draw)
		view % initial viewing angle
		waxis % initial axis
		computeH % whether to compute the energy
		plotH % whether to plot the energy at the end
		Hexpected % expected energy
		task % task for the adjoint method
		solverInfo % solver info for the current step
	end
	
	methods
		function this = Scene()
			% Default values
			redmax.Scene.clear();
			this.name = '';
			this.bodies = {};
			this.joints = {};
			this.forces = {};
			this.tEnd = 1;
			this.qInit = [];
			this.qdotInit = [];
			this.h = 1e-2;
			this.t = 0;
			this.k = 0;
			this.T0 = 0;
			this.V0 = 0;
			this.history = [];
			this.nsteps = 0;
			this.grav = [0 0 -980]';
			this.drawHz = 15;
			this.view = 3;
			this.waxis = [];
			this.computeH = true;
			this.plotH = true;
			this.Hexpected = zeros(1,2);
			this.task = [];
			this.solverInfo = [];
		end
		
		function init(this)
			if usejava('jvm')
				colormap('default'); % Restore colormap mode
			end
			njoints = length(this.joints);

			% Leaf-to-root ordering is better for reducing matrix fill
			order = this.joints{1}.getTraversalOrder();
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
			
			% Update joints
			this.joints{1}.update();
			[this.qInit,this.qdotInit] = this.joints{1}.getQ();
			
			% Initialize bodies
			nbodies = length(this.bodies);
			this.bodies = this.bodies(order); % same ordering as joints
			for i = 1 : nbodies
				this.bodies{i}.computeInertia();
				if i < nbodies
					this.bodies{i}.next = this.bodies{i+1}; %#ok<*SAGROW>
				end
			end
			
			% Initialize forces
			this.forces{1}.init();
			
			% Other initial values
			this.nsteps = ceil(this.tEnd/this.h);
			this.reset();
		end
		
		%%
		function reset(this)
			this.joints{1}.setQ(this.qInit,this.qdotInit);
			this.t = 0;
			this.k = 0;
			[this.T0,this.V0] = this.joints{1}.computeEnergies(this.grav);
			this.V0 = this.forces{1}.computeEnergy(this.V0);
			if ~isempty(this.task)
				this.task.P = 0;
			end
		end
		
		%%
		function saveHistory(this,Hl,Hu,Hp,M,f,K,D,J,dJdq)
			[q,qdot] = this.joints{1}.getQ();
			this.history(this.k).q = q;
			this.history(this.k).qdot = qdot;
			% TODO: may need to store other info (eg Euler angle chart)
			if ~isempty(this.task)
				% Used by the adjoint method
				this.history(this.k).Hl = Hl;
				this.history(this.k).Hu = Hu;
				this.history(this.k).Hp = Hp;
				this.history(this.k).M = M;
				this.history(this.k).f = f;
				this.history(this.k).K = K;
				this.history(this.k).D = D;
				this.history(this.k).J = J;
				this.history(this.k).dJdq = dJdq;
			end
			if ~isempty(this.task)
				% Used by the adjoint method
				this.task.calcStep();
			end
			if this.computeH
				[T,V] = this.joints{1}.computeEnergies(this.grav);
				V = this.forces{1}.computeEnergy(V);
				this.history(this.k).T = T;
				this.history(this.k).V = V;
				this.history(this.k).t = this.t;
			end
		end
		
		%%
		function plotEnergies(this,itype)
			if this.computeH
				T = [this.T0,this.history.T];
				V = [this.V0,this.history.V];
				t = [0,this.history.t]; %#ok<PROPLC>
				V = V - V(1);
				H = T + V;
				if this.Hexpected(itype) ~= 0
					dH = H(end) - this.Hexpected(itype);
					if abs(dH) > 1e-2
						fprintf(2,'### FAIL: %.16e ###\n',H(end));
					else
						fprintf(1,'### PASS ###\n');
					end
				end
				if this.plotH
					clf;
					plot(t,T,'.-',t,V,'.-',t,H,'.-'); %#ok<PROPLC>
					a = axis();
					a(2) = t(end); %#ok<PROPLC>
					axis(a);
					grid on;
					legend('T','V','H');
					xlabel('Time');
					ylabel('Energy');
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
			if this.drawHz > 0 && (floor(this.t*this.drawHz) > floor((this.t-this.h)*this.drawHz))
				cla;
				hold on;
				this.bodies{1}.draw();
				this.joints{1}.draw();
				this.forces{1}.draw();
				if ~isempty(this.task)
					this.task.draw();
				end
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
			[J,Jdot,dJdq,dJdotdq] = joint0.computeJacobian();
			[Mm,fm,Km,Dm] = body0.computeMassGrav(this.grav);
			[fm,Km,Dm] = body0.computeForce(fm,Km,Dm);
			[fr,Kr,Dr] = joint0.computeForce();
			
			% Inertia
			M = J'*Mm*J;
			dMdq = zeros(nr,nr,nr);
			for i = 1 : nr
				tmp = J'*Mm*dJdq(:,:,i);
				dMdq(:,:,i) = tmp' + tmp;
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
				[J_,Jdot_] = joint0.computeJacobian();
				joint0.setQ(q);
				joint0.update();
				dJdq_(:,:,i) = (J_ - J)/sqrteps;
				dJdotdq_(:,:,i) = (Jdot_ - Jdot)/sqrteps;
			end
			redmax.Scene.printError('dJ/dq',dJdq_,dJdq);
			redmax.Scene.printError('dJdot/dq',dJdotdq_,dJdotdq);
			
			% Test inertia
			dMdq_ = zeros(nr,nr,nr);
			for i = 1 : nr
				q_ = q;
				q_(i) = q_(i) + sqrteps;
				joint0.setQ(q_);
				joint0.update();
				J_ = joint0.computeJacobian();
				joint0.setQ(q);
				joint0.update();
				M_ = J_'*Mm*J_;
				dMdq_(:,:,i) = (M_ - M)/sqrteps;
			end
			redmax.Scene.printError('dM/dq',dMdq_,dMdq);
			
			% Test quadratic velocity vector
			Kqvv_ = zeros(nr,nr);
			Dqvv_ = zeros(nr,nr);
			for i = 1 : nr
				% Modify q
				q_ = q;
				q_(i) = q_(i) + sqrteps;
				joint0.setQ(q_);
				joint0.update();
				[J_,Jdot_] = joint0.computeJacobian();
				joint0.setQ(q);
				joint0.update();
				fqvv_ = -J_'*Mm*Jdot_*qdot;
				Kqvv_(:,i) = (fqvv_ - fqvv)/sqrteps;
				% Modify qdot
				qdot_ = qdot;
				qdot_(i) = qdot_(i) + sqrteps;
				joint0.setQdot(qdot_);
				joint0.update();
				[J_,Jdot_] = joint0.computeJacobian();
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
				[J_,Jdot_] = joint0.computeJacobian();
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
				[J_,Jdot_] = joint0.computeJacobian();
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
