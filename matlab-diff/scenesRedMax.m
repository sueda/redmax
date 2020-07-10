function scene = scenesRedMax(sceneID)
% Creates test scenes

BDF1 = 1;
BDF2 = 2;

scene = redmax.Scene();

% Default values
density = 1.0;

switch(sceneID)
	case -1
		scene.name = 'Simpler serial chain';
		%scene.grav = [0 0 0]';
		%scene.tspan = [0 0.1];
		scene.drawHz = 15;
		sides = [10 1 1];
		nbodies = 2;
		scene.waxis = nbodies*5*[-1 1 -1 1 -2 0];
		for i = 1 : nbodies
			scene.bodies{i} = redmax.BodyCuboid(density,sides); %#ok<*SAGROW>
			if i == 1
				scene.joints{i} = redmax.JointRevolute([],scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform(eye(4));
				scene.joints{i}.q(1) = 0;%pi/2;
				scene.joints{i}.qdot(1) = 0;%1;
			else
				scene.joints{i} = redmax.JointRevolute(scene.joints{i-1},scene.bodies{i},[0 1 0]');
				scene.joints{i}.setJointTransform([eye(3),[10 0 0]'; 0 0 0 1]);
				scene.joints{i}.q(1) = pi/4;
				scene.joints{i}.qdot(1) = 1;
			end
			scene.bodies{i}.setBodyTransform([eye(3),[5 0 0]'; 0 0 0 1]);
			% Joint stiffness and damping
			scene.joints{i}.setStiffness(1e6);
			scene.joints{i}.setDamping(1e4);
		end
	case 0
		scene.name = 'Simple serial chain';
		scene.Hexpected(BDF1) = -1.2807395786125271e+05;
		scene.Hexpected(BDF2) = 1.7541950565959269e+03;
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
		scene.Hexpected(BDF1) = -3.8520013298836180e+04;
		scene.Hexpected(BDF2) = -9.0813101725242450e+02;
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
		scene.Hexpected(BDF1) = -2.2910534778856141e+04;
		scene.Hexpected(BDF2) = -2.5064584562690561e+02;
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
	case 3
		scene.name = 'Prismatic joint';
		scene.Hexpected(BDF1) = -3.7748040716659023e+04;
		scene.Hexpected(BDF2) = -8.2848134912282694e+02;
		scene.waxis = 15*[-1.0 1.0 -0.2 0.2 -1 0.1];
		scene.bodies{1} = redmax.BodyCuboid(density,[20 1 1]);
		scene.joints{1} = redmax.JointPrismatic([],scene.bodies{1},[1 0 0]');
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform(eye(4));
		scene.bodies{2} = redmax.BodyCuboid(density,[1 1 10]);
		scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[0 1 0]');
		scene.joints{2}.setJointTransform([eye(3),[-10 0 0]'; 0 0 0 1]);
		scene.bodies{2}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.joints{2}.q = pi/2;
	case 4
		scene.name = 'Planar joint';
		scene.Hexpected(BDF1) = -4.6004266902085845e+04;
		scene.Hexpected(BDF2) = -4.9890367855813383e+02;
		scene.waxis = 15*[-1.0 1.0 -1 1 -1 0.1];
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
	case 5
		scene.name = 'Translational joint';
		scene.Hexpected(BDF1) = 3.3640234400924252e+04;
		scene.Hexpected(BDF2) = 3.3374768837480864e+04;
		scene.tspan = [0.0 2.0];
		scene.grav = [0 0 0]';
		scene.waxis = 20*[-1.5 0.5 -1.0 1.0 -0.5 0.5];
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
		scene.joints{2}.qdot(1) = -10.0;
		scene.joints{3}.qdot(1) = 10.0;
	case 6
		scene.name = 'Free2D joint';
		scene.Hexpected(BDF1) = 1.4078333333333315e+01;
		scene.Hexpected(BDF2) = 1.4583797058590157e+01;
		scene.h = 1e-1;
		scene.tspan = [0.0 10.0];
		scene.grav = [0 -1 0]';
		scene.view = 2;
		scene.waxis = 10*[-1 1 -1 1 -0.1 0.1];
		scene.bodies{1} = redmax.BodyCuboid(density,[1 1 1]);
		scene.joints{1} = redmax.JointFree2D([],scene.bodies{1});
		scene.joints{1}.q = [-10 -10 0]';
		scene.joints{1}.qdot = [2 5 1]';
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform(eye(4));
	case 7
		scene.name = 'Spherical joint';
		scene.Hexpected(BDF1) = -8.8261586699291947e+03;
		scene.Hexpected(BDF2) = 8.6566316039286758e+03;
		scene.tspan = [0.0 1.0];
		scene.h = 2e-3; % double pendulum requires small time steps
		scene.drawHz = 15;
		sides = [1 1 10];
		scene.waxis = 10*[-1 1 -1 1 -1.9 0.1];
		scene.bodies{1} = redmax.BodyCuboid(density,sides);
		scene.bodies{1}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.joints{1} = redmax.JointSpherical([],scene.bodies{1}); %#ok<*UNRCH>
		scene.joints{1}.setJointTransform(eye(4));
		scene.joints{1}.q = redmax.JointSpherical.getEulerInv(scene.joints{1}.chart,se3.aaToMat([1 0 0],pi/8));
		scene.joints{1}.qdot = [2 2 2]';
		scene.bodies{2} = redmax.BodyCuboid(density,sides);
		scene.bodies{2}.setBodyTransform([eye(3),[0 0 -5]'; 0 0 0 1]);
		scene.joints{2} = redmax.JointSpherical(scene.joints{1},scene.bodies{2});
		%scene.joints{2} = redmax.JointRevolute(scene.joints{1},scene.bodies{2},[1 0 0]');
		scene.joints{2}.setJointTransform([eye(3),[0 0 -10]'; 0 0 0 1]);
		scene.joints{2}.q(1) = pi/2;
	case 8
		scene.name = 'Loop';
		scene.waxis = [-15 15 -1 1 -20 2];
		scene.Hexpected(BDF1) = 1.2346847548979090e+03;
		scene.Hexpected(BDF2) = 4.1269153806366667e+03;
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
		scene.forces{1} = redmax.ForcePointPoint(scene.bodies{3},[0 0 -5]',scene.bodies{4},[10 0 0]');
		scene.forces{1}.setStiffness(1e8);
		scene.joints{5}.qdot(1) = 5;
	case 9
		scene.name = 'Free3D joint';
		scene.Hexpected(BDF1) = 4.3958752811887525e+00;
		scene.Hexpected(BDF2) = 4.5472012875148602e+00;
		scene.h = 5e-2;
		scene.tspan = [0.0 6.0];
		scene.grav = [0 0 -1]';
		scene.waxis = 5*[-0.2 0.2 -0.2 0.2 -1.0 1.0];
		scene.bodies{1} = redmax.BodyCuboid(density,[1 1 1]);
		scene.joints{1} = redmax.JointFree3D([],scene.bodies{1});
		scene.joints{1}.qdot = [ 0 0 3 0.2 0.4 0.6]';
		scene.joints{1}.setJointTransform(eye(4));
		scene.bodies{1}.setBodyTransform(eye(4));
end

end
