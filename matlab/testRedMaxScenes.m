function scene = testRedMaxScenes(itype,sceneID)
% Creates test scenes

global RECURS_ODE45 REDMAX_ODE45 REDMAX_EULER

scene = redmax.Scene();

density = 1.0;
vel = (itype == REDMAX_EULER); % Velocity level integrator?
switch(sceneID)
	case -1
		scene.name = 'Simpler serial chain';
		sides = [10 1 1];
		nbodies = 2;
		scene.waxis = nbodies*5*[-1 1 -1 1 -2 0];
		scene.Hexpected(RECURS_ODE45) = -5.6531026717020723;
		scene.Hexpected(REDMAX_ODE45) = -5.6531026765951538;
		scene.Hexpected(REDMAX_EULER) = -3697.4545694454454861;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform(eye(4));
				scene.joints{i}.q(1) = pi/2;
			else
				scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
				scene.joints{i}.q(1) = pi/4;
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		end
	case 0
		scene.name = 'Simple serial chain';
		sides = [10 1 1];
		nbodies = 5;
		scene.waxis = nbodies*5*[-1 1 -0.1 0.1 -2 0];
		scene.Hexpected(RECURS_ODE45) = -3.0971281943493523;
		scene.Hexpected(REDMAX_ODE45) = -3.0971281068341341;
		scene.Hexpected(REDMAX_EULER) = -5930.8171118834870867;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform(eye(4));
			else
				if mod(i,2) == 1
					scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				else
					scene.joints{i} = redmax.JointFixed(scene.joints{i-1},scene.bodies{i});
				end
				%scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
			if mod(i,2) == 1
				scene.joints{i}.q(1) = pi/4;
			else
				scene.joints{i}.q(1) = 0;
			end
		end
	case 1
		scene.name = 'Different revolute axes';
		sides = [10 1 1];
		scene.waxis = [-20 20 -20 20 -20 5];
		scene.Hexpected(RECURS_ODE45) = -1.9548841516880202;
		scene.Hexpected(REDMAX_ODE45) = -1.9548841526830074;
		scene.Hexpected(REDMAX_EULER) = -9423.2594023734018265;
		scene.bodies{1} = redmax.BodyCuboid(density,sides);
		scene.bodies{2} = redmax.BodyCuboid(density,sides);
		scene.bodies{3} = redmax.BodyCuboid(density,sides);
		scene.joints{1} = redmax.JointRevolute([],       scene.bodies{1},[0 0 1]');
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{3} = redmax.JointRevolute(scene.joints{2},scene.bodies{3},[0 0 1]');
		scene.bodies{1}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		scene.bodies{2}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		scene.bodies{3}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		scene.joints{1}.setJointTransform(eye(4));
		scene.joints{2}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
		scene.joints{3}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
		scene.joints{1}.q(1) = 0;
		scene.joints{2}.q(1) = pi/2;
		scene.joints{3}.q(1) = pi/2;
	case 2
		scene.name = 'Branching';
		%      X        in YZ plane
		%      |        
		%  X---Z---Y    joint axes shown
		%  |       |
		%  |       |
		scene.waxis = [-15 15 -15 15 -10 20];
		scene.Hexpected(RECURS_ODE45) = -3.2850447686942061;
		scene.Hexpected(REDMAX_ODE45) = -3.2850447782984702;
		scene.Hexpected(REDMAX_EULER) = -1123.9825362491046690;
		scene.bodies{1} = redmax.BodyCuboid(density,[1  1 10]);
		scene.bodies{2} = redmax.BodyCuboid(density,[1 20  1]);
		scene.bodies{3} = redmax.BodyCuboid(density,[1  1 10]);
		scene.bodies{4} = redmax.BodyCuboid(density,[1  1 10]);
		scene.joints{1} = redmax.JointRevolute([],       scene.bodies{1},[1 0 0]');
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 0 1]');
		scene.joints{3} = redmax.JointRevolute(scene.joints{2},scene.bodies{3},[1 0 0]');
		scene.joints{4} = redmax.JointRevolute(scene.joints{2},scene.bodies{4},[0 1 0]');
		scene.bodies{1}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.bodies{2}.setBodyTransform([eye(3),[0 0  0]'; 0 0 0 1]);
		scene.bodies{3}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.bodies{4}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.joints{1}.setJointTransform([eye(3),[0   0  15]'; 0 0 0 1]);
		scene.joints{2}.setJointTransform([eye(3),[0   0 -10]'; 0 0 0 1]);
		scene.joints{3}.setJointTransform([eye(3),[0 -10   0]'; 0 0 0 1]);
		scene.joints{4}.setJointTransform([eye(3),[0  10   0]'; 0 0 0 1]);
		scene.joints{1}.q(1) = 0;
		scene.joints{2}.q(1) = 0;
		scene.joints{3}.q(1) = pi/4;
		scene.joints{4}.q(1) = pi/4;
	case 3
		scene.name = 'Spherical joint';
		if itype == REDMAX_EULER
			scene.tspan = [0.0 3.0];
		else
			% Reparameterizing with ode45 takes too much effort
			scene.tspan = [0.0 0.0];
		end
		sides = [1 1 10];
		scene.waxis = 10*[-1 1 -1 1 -1.9 0.1];
		scene.Hexpected(RECURS_ODE45) = 0;
		scene.Hexpected(REDMAX_ODE45) = 0;
		scene.Hexpected(REDMAX_EULER) = 7788.8055603543098186;
		scene.bodies{1} = redmax.BodyCuboid(density,sides);
		scene.bodies{1}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.joints{1} = redmax.JointSphericalExp([],scene.bodies{1}); %#ok<*UNRCH>
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{2} = redmax.BodyCuboid(density,sides);
		scene.bodies{2}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.joints{2} = redmax.JointSphericalExp(scene.joints{1},scene.bodies{2});
		scene.joints{2}.setJointTransform([eye(3),[0 0 -10]'; 0 0 0 1]);
		scene.joints{1}.q(1) = pi/8;
		scene.joints{1}.q(3) = 0;
		scene.joints{1}.qdot(3) = 2;
		scene.joints{2}.q(1) = pi/8;
	case 4
		scene.name = 'Loop';
		scene.waxis = [-15 15 -1 1 -20 2];
		scene.Hexpected(RECURS_ODE45) = 4176.3993502426255873;
		scene.Hexpected(REDMAX_ODE45) = 4176.3993502425073530;
		scene.Hexpected(REDMAX_EULER) = 3987.2011847696289806;
		scene.bodies{1} = redmax.BodyCuboid(density,[20 1  1]);
		scene.bodies{2} = redmax.BodyCuboid(density,[ 1 1 10]);
		scene.bodies{3} = redmax.BodyCuboid(density,[ 1 1 10]);
		scene.bodies{4} = redmax.BodyCuboid(density,[20 1  1]);
		scene.bodies{5} = redmax.BodyCuboid(density,[ 1 1 10]);
		scene.joints{1} = redmax.JointFixed([],scene.bodies{1});
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{3} = redmax.JointRevolute(scene.joints{1},scene.bodies{3},[0 1 0]');
		scene.joints{4} = redmax.JointRevolute(scene.joints{2},scene.bodies{4},[0 1 0]');
		scene.joints{5} = redmax.JointRevolute(scene.joints{4},scene.bodies{5},[0 1 0]');
		scene.bodies{1}.setBodyTransform(eye(4));
		scene.bodies{2}.setBodyTransform([eye(3),[ 0 0 -5]'; 0 0 0 1]);
		scene.bodies{3}.setBodyTransform([eye(3),[ 0 0 -5]'; 0 0 0 1]);
		scene.bodies{4}.setBodyTransform([eye(3),[10 0  0]'; 0 0 0 1]);
		scene.bodies{5}.setBodyTransform([eye(3),[ 0 0 -5]'; 0 0 0 1]);
		scene.joints{1}.setJointTransform(eye(4));
		scene.joints{2}.setJointTransform([eye(3),[-10 0   0]'; 0 0 0 1]);
		scene.joints{3}.setJointTransform([eye(3),[ 10 0   0]'; 0 0 0 1]);
		scene.joints{4}.setJointTransform([eye(3),[  0 0 -10]'; 0 0 0 1]);
		scene.joints{5}.setJointTransform([eye(3),[ 10 0   0]'; 0 0 0 1]);
		scene.constraints{1} = redmax.ConstraintLoop(scene.bodies{3},scene.bodies{4});
		scene.constraints{1}.setPositions([0 0 -5]',[10 0 0]');
		scene.joints{5}.qdot(1) = 5;
	case 5
		scene.name = 'Joint torque';
		scene.tspan = [0 10];
		scene.hEuler = 5e-2;
		scene.grav = [0 0 0]';
		sides = [10 1 1];
		nbodies = 3;
		scene.waxis = [];
		scene.Hexpected(RECURS_ODE45) = 160.820781707015;
		scene.Hexpected(REDMAX_ODE45) = 160.820781710469;
		scene.Hexpected(REDMAX_EULER) = 170.5971183034905607;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform(eye(4));
			else
				scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		end
		scene.sceneFcn = @sceneFcn05;
	case 6
		scene.name = 'Joint limits';
		if itype ~= REDMAX_EULER
			scene.tspan = [0 0];
		end
		sides = [10 1 1];
		nbodies = 2;
		scene.waxis = nbodies*5*[-2.5 2.5 -0.1 0.1 -2.0 0.5];
		scene.Hexpected(RECURS_ODE45) = 0;
		scene.Hexpected(REDMAX_ODE45) = 0;
		scene.Hexpected(REDMAX_EULER) = 36957.4447830002754927;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform(eye(4));
			else
				scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
			% NOTE: ode45 seems to get stuck if all the scene.constraints become
			% active at the same time. So, it's safer to not apply joint
			% limits to the root joint.
			if i > 1
				scene.constraints{end+1} = redmax.ConstraintJointLimit(scene.joints{i});
				scene.constraints{end}.setLimits(-pi/4,pi/4);
			end
		end
	case 7
		scene.name = 'Equality constrained angles';
		scene.hEuler = 2e-2;
		sides = [10 1 1];
		nbodies = 3;
		scene.waxis = nbodies*5*[-2.5 2.5 -1 1 -2.0 0.5];
		scene.Hexpected(RECURS_ODE45) = -229.7121594171621837;
		scene.Hexpected(REDMAX_ODE45) = -229.7121593945485074;
		scene.Hexpected(REDMAX_EULER) = 42645.1541420989669859;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform(eye(4));
			else
				scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
			if i > 1
				scene.constraints{end+1} = redmax.ConstraintMultQ(scene.joints{i-1},scene.joints{i});
				scene.constraints{end}.setFactor(0.5);
			end
		end
	case 8
		scene.name = 'Equality and loop';
		scene.hEuler = 2e-2;
		scene.waxis = [-15 15 -5 5 -30 1];
		scene.Hexpected(RECURS_ODE45) = 16707.8585386558770551;
		scene.Hexpected(REDMAX_ODE45) = 16707.8585386558261234;
		scene.Hexpected(REDMAX_EULER) = 14677.4348748325592169;
		scene.bodies{1} = redmax.BodyCuboid(density,[10 1  1]);
		scene.bodies{2} = redmax.BodyCuboid(density,[ 1 1 10]);
		scene.bodies{3} = redmax.BodyCuboid(density,[ 1 1 10]);
		scene.bodies{4} = redmax.BodyCuboid(density,[10 1  1]);
		scene.bodies{5} = redmax.BodyCuboid(density,[ 1 1 10]);
		scene.bodies{6} = redmax.BodyCuboid(density,[ 1 1 10]);
		scene.bodies{7} = redmax.BodyCuboid(density,[ 1 1 10]);
		scene.joints{1} = redmax.JointFixed([],scene.bodies{1});
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{3} = redmax.JointRevolute(scene.joints{2},scene.bodies{3},[0 1 0]');
		scene.joints{4} = redmax.JointRevolute(scene.joints{3},scene.bodies{4},[0 1 0]');
		scene.joints{5} = redmax.JointRevolute(scene.joints{4},scene.bodies{5},[0 1 0]');
		scene.joints{6} = redmax.JointRevolute(scene.joints{5},scene.bodies{6},[0 1 0]');
		scene.joints{7} = redmax.JointRevolute(scene.joints{4},scene.bodies{7},[0 1 0]');
		scene.bodies{1}.setBodyTransform(eye(4));
		scene.bodies{2}.setBodyTransform([eye(3),[ 0 0 -5]'; 0 0 0 1]);
		scene.bodies{3}.setBodyTransform([eye(3),[ 0 0 -5]'; 0 0 0 1]);
		scene.bodies{4}.setBodyTransform([eye(3),[-5 0  0]'; 0 0 0 1]);
		scene.bodies{5}.setBodyTransform([eye(3),[ 0 0  5]'; 0 0 0 1]);
		scene.bodies{6}.setBodyTransform([eye(3),[ 0 0  5]'; 0 0 0 1]);
		scene.bodies{7}.setBodyTransform([eye(3),[ 0 0 -5]'; 0 0 0 1]);
		scene.joints{1}.setJointTransform(eye(4));
		scene.joints{2}.setJointTransform([eye(3),[  5 0   0]'; 0 0 0 1]);
		scene.joints{3}.setJointTransform([eye(3),[  0 0 -10]'; 0 0 0 1]);
		scene.joints{4}.setJointTransform([eye(3),[  0 0 -10]'; 0 0 0 1]);
		scene.joints{5}.setJointTransform([eye(3),[-10 0   0]'; 0 0 0 1]);
		scene.joints{6}.setJointTransform([eye(3),[  0 0  10]'; 0 0 0 1]);
		scene.joints{7}.setJointTransform([eye(3),[ -5 0   0]'; 0 0 0 1]);
		scene.constraints{1} = redmax.ConstraintLoop(scene.bodies{6},scene.bodies{1});
		scene.constraints{1}.setPositions([0 0 5]',[-5 0 0]');
		scene.constraints{2} = redmax.ConstraintMultQ(scene.joints{3},scene.joints{6});
		scene.constraints{2}.setFactor(0.5);
		scene.joints{7}.qdot = 1e1;
	case 9
		% Reduced ode45 does not match redmax ode45. Ignore for now!
		scene.name = 'Hybrid dynamics';
		scene.tspan = [0 10];
		scene.hEuler = 2e-2;
		scene.grav = [0 0 0]';
		sides = [10 1 1];
		nbodies = 3;
		scene.waxis = nbodies*10*[-1 1 -0.1 0.1 -1 1];
		scene.Hexpected(RECURS_ODE45) = 27616.6067502096;
		scene.Hexpected(REDMAX_ODE45) = 29802.0819629936;
		scene.Hexpected(REDMAX_EULER) = 199570.9300431804149412;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform(eye(4));
			else
				scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		end
		if itype ~= RECURS_ODE45
			scene.constraints{end+1} = redmax.ConstraintPrescJoint(scene.joints{1},vel);
		end
		scene.sceneFcn = @sceneFcn09;
	case 10
		scene.name = 'External world force';
		scene.grav = [0 0 0]';
		sides = [10 1 1];
		nbodies = 3;
		scene.waxis = nbodies*10*[-1 1 -0.1 0.1 -0.1 1];
		scene.Hexpected(RECURS_ODE45) = 1210.7099042740396726;
		scene.Hexpected(REDMAX_ODE45) = 1210.7099042740403547;
		scene.Hexpected(REDMAX_EULER) = 1088.3425711375120954;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform(eye(4));
			else
				scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
			scene.bodies{i}.setDamping(1e1);
		end
		scene.springs{end+1} = redmax.SpringPointDirection(scene.bodies{end},[5 0 0]');
		scene.springs{end}.setDirection([0 0 1]');
		scene.springs{end}.setStiffness(1e3);
		scene.sceneFcn = @sceneFcn10;
	case 11
		scene.name = 'Joint stiffness and damping';
		scene.tspan = [0 5];
		scene.hEuler = 2e-2;
		scene.grav = [0 0 0]';
		sides = [10 1 1];
		nbodies = 3;
		stiffness = 1e4;
		damping = 1e3;
		scene.waxis = nbodies*10*[-0.1 1 -0.1 0.1 -1 1];
		scene.Hexpected(RECURS_ODE45) = 2898.56113448227;
		scene.Hexpected(REDMAX_ODE45) = 2898.56113448227;
		scene.Hexpected(REDMAX_EULER) = 2659.7218894234238178;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform(eye(4));
			else
				scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
			scene.joints{i}.setStiffness(stiffness);
			scene.joints{i}.setDamping(damping);
		end
		scene.joints{1}.qdot(1) = 1;
	case 12
		% This scene will go crazy without damping!
		% This scene doesn't work with RHD!
		scene.name = 'Mass-springs';
		if itype == RECURS_ODE45
			scene.tspan = [0.0 0.0];
		else
			scene.tspan = [0.0 1.0];
		end
		scene.hEuler = 5e-3;
		sides = [10 1 1];
		nbodies = 2;
		stiffness = 1e5;
		scene.waxis = nbodies*10*[-0.5 1.5 -0.5 0.5 -1 1];
		scene.Hexpected(RECURS_ODE45) = 0;
		scene.Hexpected(REDMAX_ODE45) = -0.0345395920267038;
		scene.Hexpected(REDMAX_EULER) = -11740.4013565295099397;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform(eye(4));
			else
				scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		end
		m = 0.1*prod(sides)*density;
		scene.deformables{end+1} = redmax.DeformableSpring(3);
		scene.deformables{end}.setStiffness(stiffness);
		scene.deformables{end}.setMass(m);
		scene.deformables{end}.setAttachments([],[10*nbodies+10,0,10]',scene.bodies{end},[5 0 0]');
		scene.deformables{end+1} = redmax.DeformableSpring(2);
		scene.deformables{end}.setStiffness(stiffness);
		scene.deformables{end}.setMass(m);
		scene.deformables{end}.setAttachments(scene.bodies{1},[0 0 0]',scene.bodies{end},[0 0 0]');
	case 13
		% This scene doesn't work with RHD!
		scene.name = 'Maximal hybrid dynamics';
		if itype == RECURS_ODE45
			scene.tspan = [0.0 0.0];
		else
			scene.tspan = [0.0 10.0];
		end
		scene.hEuler = 5e-2;
		scene.waxis = 15*[-0.5 1.5 -0.5 0.5 -1 0.5];
		scene.Hexpected(RECURS_ODE45) = 0;
		scene.Hexpected(REDMAX_ODE45) = 18805.7787972479818563;
		scene.Hexpected(REDMAX_EULER) = -765.6565884021354123;
		scene.bodies{1} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{1} = redmax.JointRevolute([],scene.bodies{1},[0 1 0]');
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.bodies{2} = redmax.BodyCuboid(density,[10 1 1]);
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{2}.setJointTransform([eye(3),[0 0 -10]'; 0 0 0 1]);
		scene.bodies{2}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		scene.bodies{3} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{3} = redmax.JointRevolute(scene.joints{2},scene.bodies{3},[0 1 0]');
		scene.joints{3}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
		scene.bodies{3}.setBodyTransform([eye(3),[0 0 5]'; 0 0 0 1]);
		scene.bodies{4} = redmax.BodyCuboid(density,[10 1 1]);
		scene.joints{4} = redmax.JointRevolute(scene.joints{3},scene.bodies{4},[0 1 0]');
		scene.joints{4}.setJointTransform([eye(3),[0 0 10]'; 0 0 0 1]);
		scene.bodies{4}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		if itype ~= RECURS_ODE45
			scene.constraints{end+1} = redmax.ConstraintPrescBody(scene.bodies{end},[2 4 6],vel);
		end
		scene.sceneFcn = @sceneFcn13;
	case 14
		scene.name = 'Universal joint';
		sides = [1 1 10];
		nbodies = 3;
		scene.waxis = nbodies*10*[-0.7 0.7 -0.7 0.7 -1 0.1];
		scene.Hexpected(RECURS_ODE45) = -0.8577782593856682;
		scene.Hexpected(REDMAX_ODE45) = -0.8577782794236555;
		scene.Hexpected(REDMAX_EULER) = 9679.3365423127470422;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides);
			if i == 1
				scene.joints{i} = redmax.JointUniversal([],scene.bodies{i});
				scene.joints{i}.setJointTransform(eye(4));
			else
				scene.joints{i} = redmax.JointUniversal(scene.joints{i-1},scene.bodies{i});
				scene.joints{i}.setJointTransform([eye(3),[0 0 -10]'; 0 0 0 1]);
			end
			scene.bodies{i}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
			if mod(i,2) == 1
				scene.joints{i}.q(1) = pi/8;
			else
				scene.joints{i}.q(2) = pi/8;
			end
		end
	case 15
		scene.name = 'Prismatic joint';
		scene.waxis = 15*[-1.0 1.0 -0.2 0.2 -1 0.1];
		scene.Hexpected(RECURS_ODE45) = 2.5092171102578504;
		scene.Hexpected(REDMAX_ODE45) = 2.5092171060550754;
		scene.Hexpected(REDMAX_EULER) = -17427.8561972516035894;
		scene.bodies{1} = redmax.BodyCuboid(density,[22 1 1]);
		scene.joints{1} = redmax.JointFixed([],scene.bodies{1});
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform(eye(4));
		scene.bodies{2} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{2}.setJointTransform([eye(3),[-11 0 0]'; 0 0 0 1]);
		scene.bodies{2}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.bodies{3} = redmax.BodyCuboid(density,[22 1 1]);
		scene.joints{3} = redmax.JointPrismatic(scene.joints{2},scene.bodies{3},[1 0 0]');
		scene.joints{3}.setJointTransform([eye(3),[0 0 -10]'; 0 0 0 1]);
		scene.bodies{3}.setBodyTransform([eye(3),[11 0 0]'; 0 0 0 1]);
		scene.joints{3}.setGeometry([10 0.5 0.5]);
		scene.bodies{4} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{4} = redmax.JointRevolute(scene.joints{3},scene.bodies{4},[0 1 0]');
		scene.joints{4}.setJointTransform([eye(3),[22 0 0]'; 0 0 0 1]);
		scene.bodies{4}.setBodyTransform([eye(3),[0 0 5]'; 0 0 0 1]);
		scene.bodies{5} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{5} = redmax.JointRevolute(scene.joints{3},scene.bodies{5},[0 1 0]');
		scene.joints{5}.setJointTransform([eye(3),[11 0 0]'; 0 0 0 1]);
		scene.bodies{5}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.constraints{1} = redmax.ConstraintLoop(scene.bodies{4},scene.bodies{1});
		scene.constraints{1}.setPositions([0 0 5]',[11 0 0]');
		scene.joints{5}.q(1) = 3*pi/4;
	case 16
		scene.name = 'Planar joint';
		scene.waxis = 15*[-1.0 1.0 -1 1 -1 0.1];
		scene.Hexpected(RECURS_ODE45) = -5.7644270883174613;
		scene.Hexpected(REDMAX_ODE45) = -5.7644270894088550;
		scene.Hexpected(REDMAX_EULER) = 1027.3404900101377279;
		scene.bodies{1} = redmax.BodyCuboid(density,[10 10 1]);
		scene.joints{1} = redmax.JointPlanar([],scene.bodies{1});
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform(eye(4));
		scene.bodies{2} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{2}.setJointTransform([eye(3),[-5 0 0]'; 0 0 0 1]);
		scene.bodies{2}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.bodies{3} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{3} = redmax.JointRevolute(scene.joints{1},scene.bodies{3},[1 0 0]');
		scene.joints{3}.setJointTransform([eye(3),[0 -5 0]'; 0 0 0 1]);
		scene.bodies{3}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.joints{2}.q(1) = pi/2;
		scene.joints{3}.q(1) = pi/4;
	case 17
		scene.name = 'Translational joint';
		scene.grav = [0 0 0]';
		scene.waxis = 20*[-1.5 0.5 -1.0 1.0 -0.5 0.5];
		scene.Hexpected(RECURS_ODE45) = 835.418079875333;
		scene.Hexpected(REDMAX_ODE45) = 835.418079875333;
		scene.Hexpected(REDMAX_EULER) = 836.2350063173605577;
		scene.bodies{1} = redmax.BodyCuboid(density,[10 10 1]);
		scene.joints{1} = redmax.JointTranslational([],scene.bodies{1});
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform(eye(4));
		scene.bodies{2} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{2}.setJointTransform([eye(3),[-5 0 0]'; 0 0 0 1]);
		scene.bodies{2}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.bodies{3} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{3} = redmax.JointRevolute(scene.joints{1},scene.bodies{3},[1 0 0]');
		scene.joints{3}.setJointTransform([eye(3),[0 -5 0]'; 0 0 0 1]);
		scene.bodies{3}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.joints{1}.qdot(1) = 0.0;
		scene.joints{1}.qdot(3) = 0.0;
		scene.joints{2}.qdot(1) = 2.0;
		scene.joints{3}.qdot(1) = 1.0;
	case 18
		scene.name = 'Free joint';
		scene.tspan = [0.0 7.0];
		scene.grav = [0 0 -1]';
		scene.waxis = 5*[-0.2 0.2 -0.2 0.2 -1.0 1.0];
		scene.Hexpected(RECURS_ODE45) = 4.5466342688068826;
		scene.Hexpected(REDMAX_ODE45) = 4.5466342688068924;
		scene.Hexpected(REDMAX_EULER) = 4.5116666666668817;
		scene.bodies{1} = redmax.BodyCuboid(density,[1 1 1]);
		freeType = 2;
		if freeType == 1
			scene.joints{1} = reduced.JointFree([],scene.bodies{1});
			scene.joints{1}.qdot = [0.2 0.4 0.6 0 0 3]';
		elseif freeType == 2
			% This is the best
			scene.joints{1} = reduced.JointFree3D([],scene.bodies{1});
			scene.joints{1}.qdot = [ 0 0 3 0.2 0.4 0.6]';
		elseif freeType == 3
			jointT = reduced.JointTranslational([],scene.bodies{1});
			jointS = reduced.JointSphericalExp([],scene.bodies{1});
			scene.joints{1} = reduced.JointComposite([],scene.bodies{1},jointT,jointS);
			scene.joints{1}.qdot = [0 0 3 0.2 0.4 0.6]';
		end
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform(eye(4));
	case 19
		scene.name = 'Composite joint';
		scene.hEuler = 2e-2;
		scene.waxis = 5*[-0.1 2.1 -2 2 -2 2];
		scene.Hexpected(RECURS_ODE45) = -8.7962825142149086;
		scene.Hexpected(REDMAX_ODE45) = -8.7962825142917609;
		scene.Hexpected(REDMAX_EULER) = 918.5086593280602756;
		scene.bodies{1} = redmax.BodyCuboid(density,[1 1 10]);
		revolute  = redmax.JointRevolute( [],scene.bodies{1},[1 0 0]');
		prismatic = redmax.JointPrismatic([],scene.bodies{1},[1 0 0]');
		scene.joints{1} = redmax.JointComposite([],scene.bodies{1},revolute,prismatic);
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform([eye(3),[0 0 5]'; 0 0 0 1]);
		scene.joints{1}.q(1) = 0.1;
		scene.joints{1}.qdot(2) = 1.0;
	case 20
		% This scene doesn't work with RHD!
		scene.name = 'Reduced/maximal hybrid dynamics';
		if itype == RECURS_ODE45
			scene.tspan = [0.0 0.0];
		else
			scene.tspan = [0.0 10.0];
		end
		scene.hEuler = 5e-2;
		scene.waxis = [];%15*[-0.5 1.5 -0.5 0.5 -1 0.5];
		scene.Hexpected(RECURS_ODE45) = 0;
		scene.Hexpected(REDMAX_ODE45) = 72822.5867580034246203;
		scene.Hexpected(REDMAX_EULER) = 50368.3587015155280824;
		scene.bodies{1} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{1} = redmax.JointRevolute([],scene.bodies{1},[0 1 0]');
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.bodies{2} = redmax.BodyCuboid(density,[10 1 1]);
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{2}.setJointTransform([eye(3),[0 0 -10]'; 0 0 0 1]);
		scene.bodies{2}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		scene.bodies{3} = redmax.BodyCuboid(density,[10 1 1]);
		scene.joints{3} = redmax.JointRevolute(scene.joints{2},scene.bodies{3},[0 1 0]');
		scene.joints{3}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
		scene.bodies{3}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		scene.bodies{4} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{4} = redmax.JointRevolute(scene.joints{3},scene.bodies{4},[0 1 0]');
		scene.joints{4}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
		scene.bodies{4}.setBodyTransform([eye(3),[0 0 5]'; 0 0 0 1]);
		scene.bodies{5} = redmax.BodyCuboid(density,[10 1 1]);
		scene.joints{5} = redmax.JointRevolute(scene.joints{4},scene.bodies{5},[0 1 0]');
		scene.joints{5}.setJointTransform([eye(3),[0 0 10]'; 0 0 0 1]);
		scene.bodies{5}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		if itype ~= RECURS_ODE45
			scene.constraints{end+1} = redmax.ConstraintPrescBody(scene.bodies{end},[2 4 6],vel);
			scene.constraints{end+1} = redmax.ConstraintPrescJoint(scene.joints{3},vel);
		end
		scene.sceneFcn = @sceneFcn20;
	case 21
		scene.name = 'Spline curve joint';
		scene.hEuler = 5e-3;
		scene.waxis = 15*[-1.5 1.5 -0.5 0.5 -3 0.0];
		scene.Hexpected(RECURS_ODE45) = -18.5261468157405034;
		scene.Hexpected(REDMAX_ODE45) = -18.5261468464450445;
		scene.Hexpected(REDMAX_EULER) = -30627.8479814097263443;
		scene.bodies{1} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{1} = redmax.JointRevolute([],scene.bodies{1},[0 1 0]');
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.bodies{2} = redmax.BodyCuboid(density,[10 1 1]);
		scene.joints{2} = redmax.JointSplineCurve(scene.joints{1},scene.bodies{2});
		scene.joints{2}.setJointTransform([eye(3),[0 0 -10]'; 0 0 0 1]);
		scene.joints{2}.addControlFrame([se3.aaToMat([0 1 0],pi),   [-10 0  0]'; 0 0 0 1]);
		scene.joints{2}.addControlFrame([se3.aaToMat([0 1 0],pi/2), [  0 0 -2]'; 0 0 0 1]);
		scene.joints{2}.addControlFrame([se3.aaToMat([0 1 0],0),    [ 10 0  0]'; 0 0 0 1]);
		scene.joints{2}.addControlFrame([se3.aaToMat([0 1 0],-pi/2),[  0 0  2]'; 0 0 0 1]);
		scene.bodies{2}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		scene.bodies{3} = redmax.BodyCuboid(density,[10 1 1]);
		scene.joints{3} = redmax.JointRevolute(scene.joints{2},scene.bodies{3},[0 1 0]');
		scene.joints{3}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
		scene.bodies{3}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		scene.joints{2}.q(1) = 0;
		scene.joints{3}.q(1) = 15*pi/16;
	case 22
		scene.name = 'Spline surface joint';
		scene.waxis = 10*[-1.5 1.5 -1.5 1.5 -3 0];
		scene.Hexpected(RECURS_ODE45) = -1.4604474130101153;
		scene.Hexpected(REDMAX_ODE45) = -1.4604474127263529;
		scene.Hexpected(REDMAX_EULER) = 2154.9740571399888722;
		scene.bodies{1} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{1} = redmax.JointRevolute([],scene.bodies{1},[1 0 0]');
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.bodies{2} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{2} = redmax.JointSplineSurface(scene.joints{1},scene.bodies{2});
		scene.joints{2}.setJointTransform([eye(3),[0 0 -14]'; 0 0 0 1]);
		t0 = 15;
		r0 = pi/4;
		for i1 = 1 : 4
			s1 = (i1-1)/(4-1);
			x = (1-s1)*(-t0) + s1*t0;
			a = (1-s1)*(-r0) + s1*r0;
			for i2 = 1 : 4
				s2 = (i2-1)/(4-1);
				y = (1-s2)*(-t0) + s2*t0;
				z = 0.05*(x*x + y*y);
				b = (1-s1)*(-r0) + s1*r0;
				c = 0;
				scene.joints{2}.setControlFrame(i1,i2,[x y z a b c]');
			end
		end
		scene.bodies{2}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.bodies{3} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{3} = redmax.JointRevolute(scene.joints{2},scene.bodies{3},[0 1 0]');
		scene.joints{3}.setJointTransform([eye(3),[0 0 -10]'; 0 0 0 1]);
		scene.bodies{3}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.joints{1}.q(1) = pi/8;
		scene.joints{2}.q(1) = 0.5;
		scene.joints{2}.q(2) = 0.5;
		scene.joints{3}.q(1) = pi/4;
	case 23
		scene.name = 'Point-to-point spring';
		sides = [10 1 1];
		nbodies = 4;
		scene.waxis = nbodies*10*[-1 1 -0.1 0.1 -1 0.1];
		scene.Hexpected(RECURS_ODE45) = -0.2671194855411159;
		scene.Hexpected(REDMAX_ODE45) = -0.2671194856266084;
		scene.Hexpected(REDMAX_EULER) = 2125.1442936080966319;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform(eye(4));
				scene.joints{i}.q(1) = pi/2;
			else
				scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
				scene.joints{i}.q(1) = pi/16;
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
			if i > 1
				scene.springs{i-1} = redmax.SpringPointPoint(scene.bodies{i-1},[-1 0 0]',scene.bodies{i},[5 0 0]');
				scene.springs{i-1}.setStiffness(1e2);
			end
		end
	case 24
		scene.name = 'Spring damper';
		if itype == RECURS_ODE45 || itype == REDMAX_ODE45
			scene.tspan = [0.0 0.0];
		end
		sides = [10 1 1];
		scene.waxis = 3*10*[-0.1 1 -0.1 0.1 -1 0.1];
		scene.Hexpected(RECURS_ODE45) = 0;
		scene.Hexpected(REDMAX_ODE45) = 0;
		scene.Hexpected(REDMAX_EULER) = -18398.2926338097677217;
		testPCG = false;
		scene.bodies{1} = redmax.BodyCuboid(density,sides);
		if testPCG
			scene.joints{1} = redmax.JointRevolute([],scene.bodies{1},[0 1 0]');
			scene.joints{1}.setStiffness(1e9);
			scene.joints{1}.setDamping(1e3);
		else
			scene.joints{1} = redmax.JointFixed([],scene.bodies{1});
		end
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		scene.bodies{2} = redmax.BodyCuboid(density,sides);
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{2}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
		scene.joints{2}.q(1) = pi/2;
		scene.bodies{2}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		scene.bodies{3} = redmax.BodyCuboid(density,sides);
		scene.joints{3} = redmax.JointRevolute(scene.joints{2},scene.bodies{3},[0 1 0]');
		scene.joints{3}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
		scene.joints{3}.q(1) = -pi/2;
		scene.bodies{3}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		scene.springs{1} = redmax.SpringDamper(scene.bodies{1},[-2 0 -0.5]',scene.bodies{2},[1 0.5 -0.5]');
		scene.springs{1}.setStiffness(1e6);
		scene.springs{1}.setDamping(1e3);
		if testPCG
			scene.springs{2} = redmax.SpringDamper(scene.bodies{1},[2 0 -0.5]',scene.bodies{3},[1 -0.5 0.5]');
			scene.springs{2}.setStiffness(1e6);
			scene.springs{2}.setDamping(1e3);
		end
	case 25
		scene.name = 'Composite body';
		sides = [1 1 10];
		scene.waxis = 10*[-1 1 -0.1 0.1 -1.5 0.1];
		scene.Hexpected(RECURS_ODE45) = -11.2086902929768257;
		scene.Hexpected(REDMAX_ODE45) = -11.2086902930313954;
		scene.Hexpected(REDMAX_EULER) = 1261.6057602036726166;
		if false
			% With a fixed joint
			scene.bodies{1} = redmax.BodyCuboid(density,sides);
			scene.joints{1} = redmax.JointRevolute([],scene.bodies{1},[0 1 0]');
			scene.joints{1}.setJointTransform(eye(4));
			scene.bodies{1}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
			scene.bodies{2} = redmax.BodyCylinder(density,1.0,10.0);
			scene.joints{2} = redmax.JointFixed(scene.joints{1},scene.bodies{2});
			scene.joints{2}.setJointTransform([eye(3),[0 0 -10]'; 0 0 0 1]);
			scene.bodies{2}.setBodyTransform([se3.aaToMat([0 1 0],pi/2),[5 0 0]'; 0 0 0 1]);
		else
			% With a composite body
			body1a = redmax.BodyCuboid(density,sides);
			body1b = redmax.BodyCylinder(density,1.0,10.0);
			scene.bodies{1} = redmax.BodyComposite(density);
			scene.joints{1} = redmax.JointRevolute([],scene.bodies{1},[0 1 0]');
			scene.joints{1}.setJointTransform(eye(4));
			scene.bodies{1}.addBody(body1a,[eye(3),[0 0 -5]'; 0 0 0 1]);
			scene.bodies{1}.addBody(body1b,[se3.aaToMat([0 1 0],pi/2),[5 0 -10]'; 0 0 0 1]);
			E = scene.bodies{1}.computeInertiaFrame();
			scene.bodies{1}.setBodyTransform(E);
		end
	case 26
		scene.name = 'Obj body';
		scene.tspan = [0 1];
		scene.waxis = 4*[-1 1 -0.5 0.5 -1.9 0.1];
		scene.Hexpected(RECURS_ODE45) = -0.0441469434378234;
		scene.Hexpected(REDMAX_ODE45) = -0.0441469434412625;
		scene.Hexpected(REDMAX_EULER) = 59.8820887155682158;
		E0 = [eye(3),[0.5 0 -1.5]'; 0 0 0 1]; % body transform
		if false
			% With BodyCuboid
			scene.bodies{1} = redmax.BodyCuboid(density,[1 2 3]);
			scene.joints{1} = redmax.JointRevolute([],scene.bodies{1},[0 1 0]');
			scene.joints{1}.setJointTransform(eye(4));
			scene.bodies{1}.setBodyTransform(E0);
			scene.bodies{2} = redmax.BodyCuboid(density,[1 2 3]);
			scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
			scene.joints{2}.setJointTransform([eye(3),[0 0 -3]'; 0 0 0 1]);
			scene.bodies{2}.setBodyTransform(E0);
		else
			% With BodyMeshObj
			scene.bodies{1} = redmax.BodyMeshObj(density,'cuboid.obj');
			scene.joints{1} = redmax.JointRevolute([],scene.bodies{1},[0 1 0]');
			scene.joints{1}.setJointTransform(eye(4));
			scene.bodies{1}.setBodyTransform(E0 * scene.bodies{1}.E_oi);
			scene.bodies{2} = redmax.BodyMeshObj(density,'cuboid.obj');
			scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
			scene.joints{2}.setJointTransform([eye(3),[0 0 -3]'; 0 0 0 1]);
			scene.bodies{2}.setBodyTransform(E0 * scene.bodies{2}.E_oi);
		end
	case 27
		scene.name = 'Internal friction revolute';
		if itype == REDMAX_EULER
			scene.tspan = [0 1];
		else
			scene.tspan = [0 0];
		end
		scene.fric = true;
		sides = [10 1 1];
		nbodies = 2;
		angle = pi/4;
		scene.waxis = nbodies*5*[-1 1 -1 1 -2 0];
		scene.Hexpected(RECURS_ODE45) = 0;
		scene.Hexpected(REDMAX_ODE45) = 0;
		scene.Hexpected(REDMAX_EULER) = -137371.1285153437056579;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 0 1]');
				scene.joints{i}.setJointTransform([se3.aaToMat([1 0 0]',angle),[0 0 0]'; 0 0 0 1]);
				%scene.joints{i}.q = -pi/2;
			else
				scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 0 1]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
			scene.joints{i}.setGeometry(0.5,1);
		end
	case 28
		scene.name = 'Internal friction spherical';
		if itype == REDMAX_EULER
			scene.tspan = [0 1];
		else
			scene.tspan = [0 0];
		end
		scene.fric = true;
		scene.mu(1) = 5.0;
		sides = [10 1 1];
		nbodies = 2;
		angle = 0;
		scene.waxis = nbodies*5*[-1 1 -1 1 -2 0];
		scene.Hexpected(RECURS_ODE45) = 0;
		scene.Hexpected(REDMAX_ODE45) = 0;
		scene.Hexpected(REDMAX_EULER) = -184565.9459125697612762;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointSphericalExp([],scene.bodies{i});
				scene.joints{i}.setJointTransform([se3.aaToMat([1 0 0]',angle),[0 0 0]'; 0 0 0 1]);
			else
				scene.joints{i} = redmax.JointSphericalExp(scene.joints{i-1},scene.bodies{i});
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
				%scene.joints{i}.q = [0 0 pi/2]';
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
			scene.joints{i}.setGeometry(1.0);
		end
	case 29
		scene.name = 'Internal friction prismatic';
		if itype == REDMAX_EULER
			scene.tspan = [0 1.0];
		else
			scene.tspan = [0 0];
		end
		scene.fric = true;
		scene.mu(1) = 0.8;
		sides = [10 1 1];
		nbodies = 2;
		angle = pi/3;
		scene.waxis = nbodies*5*[-1 1 -1 1 -2 0];
		scene.Hexpected(RECURS_ODE45) = 0;
		scene.Hexpected(REDMAX_ODE45) = 0;
		scene.Hexpected(REDMAX_EULER) = -256391.5065969563729595;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointPrismatic([],scene.bodies{i},[1 0 0]');
				scene.joints{i}.setJointTransform([se3.aaToMat([0 1 0]',angle),[0 0 0]'; 0 0 0 1]);
			else
				scene.joints{i} = redmax.JointPrismatic(scene.joints{i-1},scene.bodies{i},[1 0 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
			scene.joints{i}.setGeometry([5,0.5,5]);
		end
	case 30
		% https://en.wikipedia.org/wiki/File:Grashof_Type_I_Four-Bar_Kinematic_Inversions.gif
		scene.name = 'Internal friction 4-bar linkage';
		if itype == REDMAX_EULER
			scene.tspan = [0 1];
			scene.hEuler = 5e-3;
			scene.fric = true;
			scene.mu(1) = 0.3;
			scene.baumgarte(3) = 1/scene.hEuler;
			scene.drawHz = 30;
		else
			scene.tspan = [0 1.5];
			scene.fric = false;
		end
		scene.waxis = [-10 6 -1 1 -5 15];
		scene.Hexpected(RECURS_ODE45) = -2.1968461897658926;
		scene.Hexpected(REDMAX_ODE45) = -2.1968459311251536;
		scene.Hexpected(REDMAX_EULER) = -14581.1508526040543074;
		scene.bodies{1} = redmax.BodyCuboid(density,[10 0.5 0.5]);
		scene.joints{1} = redmax.JointFixed([],scene.bodies{1});
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{2} = redmax.BodyCuboid(density,[4 0.5 0.5]);
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{2}.setJointTransform([eye(3),[-5 0 0]'; 0 0 0 1]);
		scene.joints{2}.setGeometry(0.5,0.5);
		scene.bodies{2}.setBodyTransform([eye(3),[2 0 0]'; 0 0 0 1]);
		% https://en.wikipedia.org/wiki/Altitude_(triangle)#Altitude_in_terms_of_the_sides
		a = 6;
		b = 12;
		c = 14;
		s = 0.5*(a + b + c);
		z = 2*sqrt(s*(s-a)*(s-b)*(s-c))/a;
		x = sqrt(14*14 - z*z);
		theta = atan2(z,x);
		scene.bodies{3} = redmax.BodyCuboid(density,[14 0.5 0.5]);
		scene.joints{3} = redmax.JointRevolute(scene.joints{2},scene.bodies{3},[0 1 0]');
		scene.joints{3}.setJointTransform([eye(3),[4 0 0]'; 0 0 0 1]);
		scene.joints{3}.setGeometry(0.5,0.5);
		scene.bodies{3}.setBodyTransform([se3.aaToMat([0 1 0],-theta),[0.5*x,0,0.5*z]'; 0 0 0 1]);
		scene.bodies{4} = redmax.BodyCuboid(density,[12 0.5 0.5]);
		scene.joints{4} = redmax.JointRevolute(scene.joints{3},scene.bodies{4},[0 1 0]');
		scene.joints{4}.setJointTransform([eye(3),[x,0,z]'; 0 0 0 1]);
		scene.joints{4}.setGeometry(0.5,0.5);
		x = x - 6;
		theta = atan2(z,x);
		scene.bodies{4}.setBodyTransform([se3.aaToMat([0 1 0],-theta),[-0.5*x,0,-0.5*z]'; 0 0 0 1]);
		scene.constraints{1} = redmax.ConstraintLoop(scene.bodies{4},scene.bodies{1});
		scene.constraints{1}.setPositions([-6 0 0]',[5 0 0]');
		scene.constraints{1}.setGeometry(0.5,0.5);
	case 31
		scene.name = 'External friction';
		if itype == RECURS_ODE45 || itype == REDMAX_ODE45
			scene.tspan = [0.0 0.0];
		else
			scene.tspan = [0 2];
		end
		sides = [10 1 1];
		nbodies = 2;
		scene.waxis = nbodies*5*[-1 1 -0.1 0.1 -2 0];
		scene.fric = true;
		scene.mu = [0.1 0.2];
		scene.SPiterMax = 100;
		scene.SPconv = 1e-3;
		scene.SPathresh = 1e-10;
		scene.SPreg = 1e-6;
		scene.baumgarte(3) = 0.1/scene.hEuler;
		scene.Hexpected(RECURS_ODE45) = 0;
		scene.Hexpected(REDMAX_ODE45) = 0;
		scene.Hexpected(REDMAX_EULER) = -90558.1346001959173009;
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides);
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform(eye(4));
				scene.joints{i}.q(1) = pi/4;
			else
				scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
				scene.joints{i}.q(1) = -pi/4;
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		end
		% Add the sphere at the end
		scene.bodies{end+1} = redmax.BodySphere(density,sides(2));
		scene.joints{end+1} = redmax.JointFixed(scene.joints{end},scene.bodies{end});
		scene.joints{end}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
		scene.bodies{end}.setBodyTransform(eye(4));
		% Add floor collision
		scene.constraints{end+1} = redmax.ConstraintFloor(scene.bodies{end});
		scene.constraints{end}.setTransform([eye(3),[0 0 -15]'; 0 0 0 1]);
	case 32
		% https://en.wikipedia.org/wiki/File:Grashof_Type_I_Four-Bar_Kinematic_Inversions.gif
		scene.name = 'Prescribed joint via maximal constraint';
		if itype == REDMAX_EULER
			scene.tspan = [0 1];
			scene.hEuler = 5e-3;
			scene.fric = false;
			scene.mu = [0.1 0.1];
			scene.baumgarte(3) = 0.1/scene.hEuler;
			scene.drawHz = 30;
		else
			scene.tspan = [0 0];
		end
		scene.waxis = [-8 8 -1 1 -5 15];
		scene.Hexpected(RECURS_ODE45) = 0;
		scene.Hexpected(REDMAX_ODE45) = 0;
		scene.Hexpected(REDMAX_EULER) = 4641.9162041538456833;
		scene.bodies{1} = redmax.BodyCuboid(density,[10 0.5 0.5]);
		scene.joints{1} = redmax.JointRevolute([],scene.bodies{1},[0 1 0]');
		scene.joints{1}.setJointTransform([eye(3),[0 0 10]'; 0 0 0 1]);
		scene.joints{1}.q = pi;
		scene.bodies{2} = redmax.BodyCuboid(density,[4 0.5 0.5]);
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{2}.setJointTransform([eye(3),[-5 0 0]'; 0 0 0 1]);
		scene.joints{2}.setGeometry(0.5,0.5);
		scene.bodies{2}.setBodyTransform([eye(3),[2 0 0]'; 0 0 0 1]);
		% https://en.wikipedia.org/wiki/Altitude_(triangle)#Altitude_in_terms_of_the_sides
		a = 6;
		b = 12;
		c = 14;
		s = 0.5*(a + b + c);
		z = 2*sqrt(s*(s-a)*(s-b)*(s-c))/a;
		x = sqrt(14*14 - z*z);
		theta = atan2(z,x);
		scene.bodies{3} = redmax.BodyCuboid(density,[14 0.5 0.5]);
		scene.joints{3} = redmax.JointRevolute(scene.joints{2},scene.bodies{3},[0 1 0]');
		scene.joints{3}.setJointTransform([eye(3),[4 0 0]'; 0 0 0 1]);
		scene.joints{3}.setGeometry(0.5,0.5);
		scene.bodies{3}.setBodyTransform([se3.aaToMat([0 1 0],-theta),[0.5*x,0,0.5*z]'; 0 0 0 1]);
		scene.bodies{4} = redmax.BodyCuboid(density,[12 0.5 0.5]);
		scene.joints{4} = redmax.JointRevolute(scene.joints{3},scene.bodies{4},[0 1 0]');
		scene.joints{4}.setJointTransform([eye(3),[x,0,z]'; 0 0 0 1]);
		scene.joints{4}.setGeometry(0.5,0.5);
		x = x - 6;
		theta = atan2(z,x);
		scene.bodies{4}.setBodyTransform([se3.aaToMat([0 1 0],-theta),[-0.5*x,0,-0.5*z]'; 0 0 0 1]);
		scene.constraints{1} = redmax.ConstraintLoop(scene.bodies{4},scene.bodies{1});
		scene.constraints{1}.setPositions([-6 0 0]',[5 0 0]');
		scene.constraints{1}.setGeometry(0.5,0.5);
		% Prescribe driver
		scene.constraints{end+1} = redmax.ConstraintPrescJointM(scene.joints{2},vel);
		scene.sceneFcn = @sceneFcn32;
	case 33
		% https://en.wikipedia.org/wiki/File:Grashof_Type_I_Four-Bar_Kinematic_Inversions.gif
		scene.name = 'External friction 4-bar linkage';
		if itype == REDMAX_EULER
			scene.tspan = [0 1];
			scene.hEuler = 5e-3;
			scene.fric = true;
			scene.mu = [0.8 0.8];
			scene.baumgarte(3) = 0.1/scene.hEuler;
			scene.drawHz = 30;
		else
			scene.tspan = [0 0];
		end
		scene.waxis = [-10 6 -1 1 -5 15];
		scene.Hexpected(RECURS_ODE45) = 0;
		scene.Hexpected(REDMAX_ODE45) = 0;
		scene.Hexpected(REDMAX_EULER) = 19598.8605086512579874;
		scene.bodies{1} = redmax.BodyCuboid(density,[10 0.5 0.5]);
		scene.joints{1} = redmax.JointFree([],scene.bodies{1});
		scene.joints{1}.setJointTransform([eye(3),[0 0 1]'; 0 0 0 1]);
		scene.bodies{2} = redmax.BodyCuboid(density,[4 0.5 0.5]);
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{2}.setJointTransform([eye(3),[-5 0 0]'; 0 0 0 1]);
		scene.joints{2}.setGeometry(0.5,0.5);
		scene.bodies{2}.setBodyTransform([eye(3),[2 0 0]'; 0 0 0 1]);
		% https://en.wikipedia.org/wiki/Altitude_(triangle)#Altitude_in_terms_of_the_sides
		a = 6;
		b = 12;
		c = 14;
		s = 0.5*(a + b + c);
		z = 2*sqrt(s*(s-a)*(s-b)*(s-c))/a;
		x = sqrt(14*14 - z*z);
		theta = atan2(z,x);
		scene.bodies{3} = redmax.BodyCuboid(density,[14 0.5 0.5]);
		scene.joints{3} = redmax.JointRevolute(scene.joints{2},scene.bodies{3},[0 1 0]');
		scene.joints{3}.setJointTransform([eye(3),[4 0 0]'; 0 0 0 1]);
		scene.joints{3}.setGeometry(0.5,0.5);
		scene.bodies{3}.setBodyTransform([se3.aaToMat([0 1 0],-theta),[0.5*x,0,0.5*z]'; 0 0 0 1]);
		scene.bodies{4} = redmax.BodyCuboid(density,[12 0.5 0.5]);
		scene.joints{4} = redmax.JointRevolute(scene.joints{3},scene.bodies{4},[0 1 0]');
		scene.joints{4}.setJointTransform([eye(3),[x,0,z]'; 0 0 0 1]);
		scene.joints{4}.setGeometry(0.5,0.5);
		x = x - 6;
		theta = atan2(z,x);
		scene.bodies{4}.setBodyTransform([se3.aaToMat([0 1 0],-theta),[-0.5*x,0,-0.5*z]'; 0 0 0 1]);
		scene.constraints{1} = redmax.ConstraintLoop(scene.bodies{4},scene.bodies{1});
		scene.constraints{1}.setPositions([-6 0 0]',[5 0 0]');
		scene.constraints{1}.setGeometry(0.5,0.5);		
		% Front sphere
		scene.bodies{end+1} = redmax.BodySphere(density,1.0);
		scene.joints{end+1} = redmax.JointFixed(scene.joints{1},scene.bodies{end});
		scene.bodies{end}.setBodyTransform([eye(3),[-5 0 0]'; 0 0 0 1]);
		scene.constraints{end+1} = redmax.ConstraintFloor(scene.bodies{end});
		% Back sphere
		scene.bodies{end+1} = redmax.BodySphere(density,1.0);
		scene.joints{end+1} = redmax.JointFixed(scene.joints{1},scene.bodies{end});
		scene.bodies{end}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
		scene.constraints{end+1} = redmax.ConstraintFloor(scene.bodies{end});
		% Middle cylinder
		scene.bodies{end+1} = redmax.BodyCylinder(density,0.5,5.0);
		scene.joints{end+1} = redmax.JointFixed(scene.joints{3},scene.bodies{end});
		scene.bodies{end}.setBodyTransform([se3.aaToMat([1 0 0],pi/2),[0 0 0]'; 0 0 0 1]);
		% Middle right sphere
		scene.bodies{end+1} = redmax.BodySphere(density,1.0);
		scene.joints{end+1} = redmax.JointFixed(scene.joints{end},scene.bodies{end});
		scene.bodies{end}.setBodyTransform([eye(3),[0 2.5 0]'; 0 0 0 1]);
		scene.constraints{end+1} = redmax.ConstraintFloor(scene.bodies{end});
		% Middle left sphere
		scene.bodies{end+1} = redmax.BodySphere(density,1.0);
		scene.joints{end+1} = redmax.JointFixed(scene.joints{end-1},scene.bodies{end});
		scene.bodies{end}.setBodyTransform([eye(3),[0 -2.5 0]'; 0 0 0 1]);
		scene.constraints{end+1} = redmax.ConstraintFloor(scene.bodies{end});
		% Prescribed driver
		scene.constraints{end+1} = redmax.ConstraintPrescJointM(scene.joints{2},vel);
		scene.sceneFcn = @sceneFcn33;
	case 34
		scene.name = 'Gears';
		% cm g s
		scene.tspan = [0 1];
		scene.hEuler = 1e-2;
		scene.drawHz = 100;
		rng(0);
		scene.waxis = [-7 7 -2 2 -7 3];
		scene.Hexpected(RECURS_ODE45) = -0.1839463800694148;
		scene.Hexpected(REDMAX_ODE45) = -0.1839463800738486;
		scene.Hexpected(REDMAX_EULER) = -39.5338848225347874;
		E0 = [se3.aaToMat([1 0 0],pi/2),[0 0 0]'; 0 0 0 1]; % body transform for gear
		% Main bar
		scene.bodies{1} = redmax.BodyCuboid(density,[1 1 6]);
		scene.joints{1} = redmax.JointFixed([],scene.bodies{1});
		scene.joints{1}.setJointTransform([se3.aaToMat([0 0 1],pi),[0 0 0]';0 0 0 1]);
		scene.bodies{1}.setBodyTransform(eye(4));
		scene.joints{1}.q = pi/2;
		% Top axle
		scene.bodies{2} = redmax.BodyCylinder(density,0.2,3.5);
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{2}.setGeometry(0.1,0.1);
		scene.joints{2}.setJointTransform([eye(3),[0 -0.25 2]';0 0 0 1]);
		scene.bodies{2}.setBodyTransform([se3.aaToMat([1 0 0],pi/2),[0 0 0]';0 0 0 1]);
		% Top gear
		scene.bodies{3} = redmax.BodyMeshObj(0.1*density,'gears.obj');
		scene.joints{3} = redmax.JointFixed(scene.joints{2},scene.bodies{3});
		scene.joints{3}.setJointTransform([eye(3),[0 -0.35 0]'; 0 0 0 1]);
		scene.bodies{3}.setBodyTransform(E0 * scene.bodies{3}.E_oi);
		% Bottom axle
		scene.bodies{4} = redmax.BodyCylinder(density,0.2,3.0);
		scene.joints{4} = redmax.JointRevolute(scene.joints{1},scene.bodies{4},[0 1 0]');
		scene.joints{4}.setGeometry(0.1,0.1);
		scene.joints{4}.setJointTransform([eye(3),[0 -0.5 -2]';0 0 0 1]);
		scene.bodies{4}.setBodyTransform([se3.aaToMat([1 0 0],pi/2),[0 0 0]';0 0 0 1]);
		% Bottom gear
		scene.bodies{5} = redmax.BodyMeshObj(0.1*density,'gears.obj');
		scene.joints{5} = redmax.JointFixed(scene.joints{4},scene.bodies{5});
		scene.joints{5}.setJointTransform([eye(3),[0 -0.1 0]'; 0 0 0 1]);
		scene.bodies{5}.setBodyTransform(E0 * scene.bodies{5}.E_oi);
		% Bottom bar
		scene.bodies{6} = redmax.BodyCuboid(density,[5 0.25 0.25]);
		scene.joints{6} = redmax.JointFixed(scene.joints{4},scene.bodies{6});
		scene.joints{6}.setJointTransform([eye(3),[0 0 0]'; 0 0 0 1]);
		scene.bodies{6}.setBodyTransform([eye(3),[0 1.5 0]'; 0 0 0 1]);
		% Top bar
		scene.bodies{7} = redmax.BodyCuboid(density,[7 0.25 0.25]);
		scene.joints{7} = redmax.JointFixed(scene.joints{2},scene.bodies{7});
		scene.joints{7}.setJointTransform([eye(3),[0 0 0]'; 0 0 0 1]);
		scene.bodies{7}.setBodyTransform([eye(3),[2 1.75 0]'; 0 0 0 1]);
		% Constraint
		scene.constraints{1} = redmax.ConstraintMultQ(scene.joints{2},scene.joints{4});
		scene.constraints{1}.setFactor(-1.0);
	case 35
		scene.name = '2D free joint';
		scene.waxis = [-20 20 -20 20 -1 1];
		scene.view = 2;
		scene.tspan = [0 10];
		scene.Hexpected(RECURS_ODE45) = 167.0835245643339135;
		scene.Hexpected(REDMAX_ODE45) = 167.0835245643319240;
		scene.Hexpected(REDMAX_EULER) = 166.9232451756938644;
		scene.grav = [0 -1 0]';
		l0 = 10;
		l1 = 10;
		% Parent
		scene.bodies{1} = redmax.BodyCuboid(density,[l0 1 1]);
		if false
			% Use composite joint
			revolute = redmax.JointRevolute([],scene.bodies{1},[0 0 1]');
			planar = redmax.JointPlanar([],scene.bodies{1});
			scene.joints{1} = redmax.JointComposite([],scene.bodies{1},planar,revolute);
		else
			% Use 2D free joint
			scene.joints{1} = redmax.JointFree2D([],scene.bodies{1});
		end
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform(eye(4));
		% Child
		scene.bodies{2} = redmax.BodyCuboid(density,[l1 1 1]);
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 0 1]');
		scene.joints{2}.setJointTransform([eye(3),[l0/2 0 0]'; 0 0 0 1]);
		scene.bodies{2}.setBodyTransform([eye(3),[l1/2 0 0]'; 0 0 0 1]);
		% Initial conditions
		scene.joints{1}.q(1:2) = [0 0]';
		scene.joints{1}.qdot(1:2) = [0 0]';
		scene.joints{1}.qdot(3) = 1;
		scene.joints{2}.qdot(1) = -1;
end

end

%%
function sceneFcn05(t,scene)
if t < 3.0
	scene.joints{1}.tau = 0;
	scene.joints{2}.tau = 0;
	scene.joints{3}.tau = 1e2;
elseif t < 6.0
	scene.joints{1}.tau = 0;
	scene.joints{2}.tau = 1e2;
	scene.joints{3}.tau = -1e2;
else
	scene.joints{1}.tau = 1e2;
	scene.joints{2}.tau = -1e2;
	scene.joints{3}.tau = 0;
end
end

%%
function sceneFcn09(t,scene)
t0 = 0.0;
t1 = 7.0;
a = 7;
b = 1.5*pi;
if t < t1
	% syms a b t t0 t1
	s = 2*((t-t0)/(t1-t0)-0.5); % -1.0 < s < 1.0
	% q(t) = b/(1+exp(-a*s))
	% dq(t) = diff(q,t);
	% ddq(t) = diff(q,t,2);
	q = b/(1+exp(-a*s));
	dq = -(2*a*b*exp(a*((2*(t - t0))/(t0 - t1) + 1)))/((t0 - t1)*(exp(a*((2*t - 2*t0)/(t0 - t1) + 1)) + 1)^2);
	ddq = (8*a^2*b*exp(2*a*((2*(t - t0))/(t0 - t1) + 1)))/((t0 - t1)^2*(exp(a*((2*t - 2*t0)/(t0 - t1) + 1)) + 1)^3) - (4*a^2*b*exp(a*((2*(t - t0))/(t0 - t1) + 1)))/((t0 - t1)^2*(exp(a*((2*t - 2*t0)/(t0 - t1) + 1)) + 1)^2);
	scene.joints{1}.presc.q = q;
	scene.joints{1}.presc.qdot = dq;
	scene.joints{1}.presc.qddot = ddq;
else
	scene.joints{1}.presc.q = b;
	scene.joints{1}.presc.qdot = 0;
	scene.joints{1}.presc.qddot = 0;
end
end

%%
function sceneFcn10(t,scene) %#ok<INUSL>
E = scene.bodies{end}.E_wi;
R = E(1:3,1:3);
G = se3.Gamma([5 0 0]');
f = [0 0 1e2]';
scene.bodies{end}.wext_i = G'*R'*f;
end

%%
function sceneFcn13(t,scene)
%
% Specifying the desired world acceleration for rotation is not
% intuitive! Let's instead specify the rotation and translation separately.
%
%    vt_w <-- target linear velocity in world space
%    wt_i <-- target angular velocity in body space
%
% Both wt_i and wtdot_i are already in body space, so they do not need to
% be transformed. Since the vt_w and vtdot_w are defined in world space, we
% must transform them.
%
%    vt_i    = R'*vt_w
%
%    vtdot_i = d/dt(R'*vt_w)
%            = R'*d/dt(vt_w) + d/dt(R')*vt_w
%            = R'*vtdot_w + ([w_i]'*R')*vt_w
%            = R'*vtdot_w + [w_i]'*vt_i
%            = R'*vtdot_w - cross(w_i,vt_i)
%
E = scene.bodies{end}.E_wi;
R = E(1:3,1:3);
phi = scene.bodies{end}.phi;
if t < 2.0
	vt_w = [0 0 0]';
	wt_i = [0 -t 0]';
	vtdot_w = [0 0 0]';
	wtdot_i = [0 -1 0]';
	% At the end of this time window, w_i = [0 -2 0]'
elseif t < 4.0
	t_ = t - 4.0; % Need to start at [0 -2 0]'
	vt_w = [0 0 0]';
	wt_i = [0 t_ 0]';
	vtdot_w = [0 0 0]';
	wtdot_i = [0 1 0]';
	% At the end of this time window, w_i = [0 0 0]'
elseif t < 6
	t_ = t - 4.0;
	vt_w = [-2*t_ 0 0]'; % Need to start at [0 0 0]'
	wt_i = [0 t_ 0]'; % Need to start at [0 0 0]'
	vtdot_w = [-2 0 0]';
	wtdot_i = [0 1 0]';
	% At the end of this time window, v_w = [-4 0 0]' and w_i = [0 2 0]'
elseif t < 8
	t_ = t - 8;
	vt_w = [2*t_ 0 0]'; % Need to start at [-4 0 0]' 
	wt_i = [0 -t_ 0]'; % Need to start at [0 2 0]'
	vtdot_w = [2 0 0]';
	wtdot_i = [0 -1 0]';
else
	vt_w = [0 0 0]';
	wt_i = [0 0 0]';
	vtdot_w = [0 0 0]';
	wtdot_i = [0 0 0]';
end
vt_i = R'*vt_w;
scene.bodies{end}.presc.qdot(1:3,1) = wt_i;
scene.bodies{end}.presc.qdot(4:6,1) = vt_i;
scene.bodies{end}.presc.qddot(1:3,1) = wtdot_i;
scene.bodies{end}.presc.qddot(4:6,1) = R'*vtdot_w - cross(phi(1:3),vt_i);
end

%%
function sceneFcn20(t,scene)
% Prescribed body
E = scene.bodies{end}.E_wi;
R = E(1:3,1:3);
phi = scene.bodies{end}.phi;
if t < 2.0
	vt_w = [0 0 0]';
	wt_i = [0 -t 0]';
	vtdot_w = [0 0 0]';
	wtdot_i = [0 -1 0]';
	% At the end of this time window, w_i = [0 -2 0]'
elseif t < 4.0
	t_ = t - 4.0; % Need to start at [0 -2 0]'
	vt_w = [0 0 0]';
	wt_i = [0 t_ 0]';
	vtdot_w = [0 0 0]';
	wtdot_i = [0 1 0]';
	% At the end of this time window, w_i = [0 0 0]'
elseif t < 6
	t_ = t - 4.0;
	vt_w = [-2*t_ 0 0]'; % Need to start at [0 0 0]'
	wt_i = [0 t_ 0]'; % Need to start at [0 0 0]'
	vtdot_w = [-2 0 0]';
	wtdot_i = [0 1 0]';
	% At the end of this time window, v_w = [-4 0 0]' and w_i = [0 2 0]'
elseif t < 8
	t_ = t - 8;
	vt_w = [2*t_ 0 0]'; % Need to start at [-4 0 0]' 
	wt_i = [0 -t_ 0]'; % Need to start at [0 2 0]'
	vtdot_w = [2 0 0]';
	wtdot_i = [0 -1 0]';
else
	vt_w = [0 0 0]';
	wt_i = [0 0 0]';
	vtdot_w = [0 0 0]';
	wtdot_i = [0 0 0]';
end
vt_i = R'*vt_w;
scene.bodies{end}.presc.qdot(1:3,1) = wt_i;
scene.bodies{end}.presc.qdot(4:6,1) = vt_i;
scene.bodies{end}.presc.qddot(1:3,1) = wtdot_i;
scene.bodies{end}.presc.qddot(4:6,1) = R'*vtdot_w - cross(phi(1:3),vt_i);
% Prescribed joint
t0 = 0.0;
t1 = 10.0;
a = 7;
b = pi/2;
% syms a b t t0 t1
s = 2*((t-t0)/(t1-t0)-0.5); % -1.0 < s < 1.0
% q(t) = b/(1+exp(-a*s))
% dq(t) = diff(q,t);
% ddq(t) = diff(q,t,2);
q = b/(1+exp(-a*s));
dq = -(2*a*b*exp(a*((2*(t - t0))/(t0 - t1) + 1)))/((t0 - t1)*(exp(a*((2*t - 2*t0)/(t0 - t1) + 1)) + 1)^2);
ddq = (8*a^2*b*exp(2*a*((2*(t - t0))/(t0 - t1) + 1)))/((t0 - t1)^2*(exp(a*((2*t - 2*t0)/(t0 - t1) + 1)) + 1)^3) - (4*a^2*b*exp(a*((2*(t - t0))/(t0 - t1) + 1)))/((t0 - t1)^2*(exp(a*((2*t - 2*t0)/(t0 - t1) + 1)) + 1)^2);
scene.joints{3}.presc.q = q;
scene.joints{3}.presc.qdot = dq;
scene.joints{3}.presc.qddot = ddq;
end

%%
function sceneFcn32(t,scene)
speed = -2*(2*pi);
q = speed*t;
scene.constraints{end}.q = q;
scene.constraints{end}.qdot = speed;
scene.constraints{end}.qddot = 0;
end

%%
function sceneFcn33(t,scene)
speed = 2*(2*pi);
q = speed*t;
scene.constraints{end}.q = q;
scene.constraints{end}.qdot = speed;
scene.constraints{end}.qddot = 0;
end
