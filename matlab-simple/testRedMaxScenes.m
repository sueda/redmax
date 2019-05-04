function scene = testRedMaxScenes(sceneID)
% Creates test scenes

scene = redmax.Scene();

density = 1.0;
switch(sceneID)
	case -1
		scene.name = 'Simpler serial chain';
		sides = [10 1 1];
		nbodies = 2;
		scene.waxis = nbodies*5*[-1 1 -1 1 -2 0];
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
			% Body damping
			scene.bodies{i}.setDamping(1e0);
			% Joint stiffness and damping
			scene.joints{i}.setStiffness(1e6);
			scene.joints{i}.setDamping(1e0);
		end
	case 0
		scene.name = 'Simple serial chain';
		sides = [10 1 1];
		nbodies = 5;
		scene.waxis = nbodies*5*[-1 1 -0.1 0.1 -2 0];
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

end

end
