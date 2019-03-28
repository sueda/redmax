classdef se3
	%se3 SE(3) and se(3) functions
	
	properties (Constant)
		THRESH = 1e-9;
	end
	
	methods (Static)

		%%
		function [ Ei ] = inv( E )
			% Inverts a transformation matrix
			R = E(1:3,1:3);
			p = E(1:3,4);
			Ei = [R', -R' * p; 0 0 0 1];
		end
		
		%%
		function [ v ] = cross( v1, v2 )
			% Matlab's version was a bottleneck!
			
			% From vecmath
			v = zeros(3,1);
			v1x = v1(1);
			v1y = v1(2);
			v1z = v1(3);
			v2x = v2(1);
			v2y = v2(2);
			v2z = v2(3);
			x = v1y*v2z - v1z*v2y;
			y = v2x*v1z - v2z*v1x;
			v(3) = v1x*v2y - v1y*v2x;
			v(1) = x;
			v(2) = y;
		end

		%%
		function [ G ] = Gamma( r )
			% Gets the 3x6 Gamma matrix, for computing the point velocity
			G = [se3.brac(r(1:3))', eye(3)];
		end
		
		%%
		function A = Ad(E)
			% Gets the adjoint transform
			A = zeros(6);
			R = E(1:3,1:3);
			p = E(1:3,4);
			A(1:3,1:3) = R;
			A(4:6,4:6) = R;
			A(4:6,1:3) = se3.brac(p) * R;
		end
		
		%%
		function a = ad(phi)
			% Gets spatial cross product matrix
			a = zeros(6);
			if size(phi,1) == 6
				w = phi(1:3,1);
				v = phi(4:6,1);
				W = se3.brac(w);
			else
				W = phi(1:3,1:3);
				v = phi(1:3,4);
			end
			a(1:3,1:3) = W;
			a(4:6,1:3) = se3.brac(v);
			a(4:6,4:6) = W;
		end
		
		%%
		function dA = Addot(E,phi)
			% Gets the time derivative of the adjoint: Ad*ad
			dA = zeros(6);
			R = E(1:3,1:3);
			p = E(1:3,4);
			w = phi(1:3);
			v = phi(4:6);
			wbrac = se3.brac(w);
			vbrac = se3.brac(v);
			pbrac = se3.brac(p);
			Rwbrac = R*wbrac;
			dA(1:3,1:3) = Rwbrac;
			dA(4:6,4:6) = Rwbrac;
			dA(4:6,1:3) = R*vbrac + pbrac*Rwbrac;
		end
		
		%%
		function [ S ] = brac( x )
			% Gets S = [x], the skew symmetric matrix
			if length(x) < 6
				S = [0, -x(3), x(2); x(3), 0, -x(1); -x(2), x(1), 0];
			else
				S(1:3,1:3) = [0, -x(3), x(2); x(3), 0, -x(1); -x(2), x(1), 0];
				S(1:3,4) = [x(4); x(5); x(6)];
				S(4,:) = zeros(1,4);
			end
		end
		
		%%
		function x = unbrac(S)
			% Gets ]x[ = S, the vector corresponding to a skew symmetric matrix
			if size(S,2) < 4
				x = [S(3,2); S(1,3); S(2,1)];
			else
				x = [S(3,2); S(1,3); S(2,1); S(1,4); S(2,4); S(3,4)];
			end
		end
		
		%%
		function R = aaToMat(axis, angle)
			% Create a rotation matrix from an (axis,angle) pair
			% From vecmath
			R = eye(3);
			ax = axis(1);
			ay = axis(2);
			az = axis(3);
			mag = sqrt( ax*ax + ay*ay + az*az);
			if mag > se3.THRESH
				mag = 1.0/mag;
				ax = ax*mag;
				ay = ay*mag;
				az = az*mag;
				if abs(ax) < se3.THRESH && abs(ay) < se3.THRESH
					% Rotation about Z
					if az < 0
						angle = -angle;
					end
					sinTheta = sin(angle);
					cosTheta = cos(angle);
					R(1,1) =  cosTheta;
					R(1,2) = -sinTheta;
					R(2,1) =  sinTheta;
					R(2,2) =  cosTheta;
				elseif abs(ay) < se3.THRESH && abs(az) < se3.THRESH
					% Rotation about X
					if ax < 0
						angle = -angle;
					end
					sinTheta = sin(angle);
					cosTheta = cos(angle);
					R(2,2) =  cosTheta;
					R(2,3) = -sinTheta;
					R(3,2) =  sinTheta;
					R(3,3) =  cosTheta;
				elseif abs(az) < se3.THRESH && abs(ax) < se3.THRESH
					% Rotation about Y
					if ay < 0
						angle = -angle;
					end
					sinTheta = sin(angle);
					cosTheta = cos(angle);
					R(1,1) =  cosTheta;
					R(1,3) =  sinTheta;
					R(3,1) = -sinTheta;
					R(3,3) =  cosTheta;
				else
					% General rotation
					sinTheta = sin(angle);
					cosTheta = cos(angle);
					t = 1.0 - cosTheta;
					xz = ax * az;
					xy = ax * ay;
					yz = ay * az;
					R(1,1) = t * ax * ax + cosTheta;
					R(1,2) = t * xy - sinTheta * az;
					R(1,3) = t * xz + sinTheta * ay;
					R(2,1) = t * xy + sinTheta * az;
					R(2,2) = t * ay * ay + cosTheta;
					R(2,3) = t * yz - sinTheta * ax;
					R(3,1) = t * xz - sinTheta * ay;
					R(3,2) = t * yz + sinTheta * ax;
					R(3,3) = t * az * az + cosTheta;
				end
			end
		end
		
		%%
		function R = qToMat(q)
			% Create a rotation matrix from a quaternion, q = [w x y z]
			% https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
			qr = q(1);
			qi = q(2);
			qj = q(3);
			qk = q(4);
			qii = qi*qi; qij = qi*qj; qik = qi*qk; qir = qi*qr;
			             qjj = qj*qj; qjk = qj*qk; qjr = qj*qr;
									  qkk = qk*qk; qkr = qk*qr;
			s = 2/(qii + qjj + qkk + qr*qr);
			R(1,1) = 1 - s*(qjj + qkk);
			R(2,1) = s*(qij + qkr);
			R(3,1) = s*(qik - qjr);
			R(1,2) = s*(qij - qkr);
			R(2,2) = 1 - s*(qii + qkk);
			R(3,2) = s*(qjk + qir);
			R(1,3) = s*(qik + qjr);
			R(2,3) = s*(qjk - qir);
			R(3,3) = 1 - s*(qii + qjj);
		end
		
		%%
		function q = matToQ(R)
			% http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
			
			m00 = R(1,1); m01 = R(1,2); m02 = R(1,3);
			m10 = R(2,1); m11 = R(2,2); m12 = R(2,3);
			m20 = R(3,1); m21 = R(3,2); m22 = R(3,3);
			tr = m00 + m11 + m22;
			if tr > 0
				S = sqrt(tr+1.0) * 2; % S=4*qw
				qw = 0.25 * S;
				qx = (m21 - m12) / S;
				qy = (m02 - m20) / S;
				qz = (m10 - m01) / S;
			elseif (m00 > m11)&&(m00 > m22)
				S = sqrt(1.0 + m00 - m11 - m22) * 2; % S=4*qx
				qw = (m21 - m12) / S;
				qx = 0.25 * S;
				qy = (m01 + m10) / S;
				qz = (m02 + m20) / S;
			elseif m11 > m22
				S = sqrt(1.0 + m11 - m00 - m22) * 2; % S=4*qy
				qw = (m02 - m20) / S;
				qx = (m01 + m10) / S;
				qy = 0.25 * S;
				qz = (m12 + m21) / S;
			else
				S = sqrt(1.0 + m22 - m00 - m11) * 2; % S=4*qz
				qw = (m10 - m01) / S;
				qx = (m02 + m20) / S;
				qy = (m12 + m21) / S;
				qz = 0.25 * S;
			end
			q = [qw qx qy qz]';
		end
		
		%%
		function E = dqToMat(dq)
			% Create a transformation matrix from a dual quaternion.
			% https://www.cs.utah.edu/~ladislav/dq/dqconv.c
			q0 = dq(:,1)/norm(dq(:,1));
			t(1,1) = 2.0*(-dq(1,2)*dq(2,1) + dq(2,2)*dq(1,1) - dq(3,2)*dq(4,1) + dq(4,2)*dq(3,1));
			t(2,1) = 2.0*(-dq(1,2)*dq(3,1) + dq(2,2)*dq(4,1) + dq(3,2)*dq(1,1) - dq(4,2)*dq(2,1));
			t(3,1) = 2.0*(-dq(1,2)*dq(4,1) - dq(2,2)*dq(3,1) + dq(3,2)*dq(2,1) + dq(4,2)*dq(1,1));
			R = se3.qToMat(q0);
			E = [R,t;0 0 0 1];
		end
		
		%%
		function dq = matToDq(E)
			% Create a transformation matrix from a dual quaternion.
			% https://www.cs.utah.edu/~ladislav/dq/dqconv.c
			t = E(1:3,4);
			q0 = se3.matToQ(E(1:3,1:3));
			dq(:,1) = q0;
			dq(1,2) = -0.5*(t(1)*q0(2) + t(2)*q0(3) + t(3)*q0(4));
			dq(2,2) = 0.5*( t(1)*q0(1) + t(2)*q0(4) - t(3)*q0(3));
			dq(3,2) = 0.5*(-t(1)*q0(4) + t(2)*q0(1) + t(3)*q0(2));
			dq(4,2) = 0.5*( t(1)*q0(3) - t(2)*q0(2) + t(3)*q0(1));
		end
		
		%%
		function E = exp(phi)
			% SE(3) exponential
			if size(phi,2) > 1
				% Convert from skew symmetric matrix to vector
				phi = se3.unbrac(phi);
			end
			% Rotational part
			I = eye(3);
			w = phi(1:3);
			wlen = norm(w);
			R = I;
			if wlen > se3.THRESH
				w = w / wlen;
				% Rodrigues formula %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				wX = w(1);
				wY = w(2);
				wZ = w(3);
				c = cos(wlen);
				s = sin(wlen);
				c1 = 1 - c;
				R = [
					c + wX*wX*c1, -wZ*s + wX*wY*c1, wY*s + wX*wZ*c1;
					wZ*s + wX*wY*c1, c + wY*wY*c1, -wX*s + wY*wZ*c1;
					-wY*s + wX*wZ*c1, wX*s + wY*wZ*c1, c + wZ*wZ*c1];
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			end
			E = R;
			% Translational part
			if length(phi) == 6
				E = [E [0 0 0]'; 0 0 0 1];
				v = phi(4:6);
				if wlen > se3.THRESH
					v = v / wlen;
					A = I - R;
					cc = se3.cross(w, v);
					d = A * cc;
					wv = w' * v;
					p = (wv * wlen) * w + d;
					E(1:3,4) = p;
				else
					E(1:3,4) = v;
				end
			end
		end
		
		%%
		function phibrac = log(E)
			% SE(3) logarithm

			% https://pixhawk.org/_media/dev/know-how/jlblanco2010geometry3d_techrep.pdf
			R = E(1:3,1:3);
			cosTheta = 0.5*(trace(R)-1);
			theta = acos(cosTheta);
			if abs(theta) < se3.THRESH
				phibrac = zeros(3);
			else
				sinTheta = sin(theta);
				wBracket = theta/(2*sinTheta)*(R-R');
				phibrac = wBracket;
			end
			if size(E,2) == 4
				phibrac = [phibrac [0 0 0]'; 0 0 0 0];
				p = E(1:3,4);
				if abs(theta) < se3.THRESH
					v = p;
				else
					V = eye(3) + (1-cosTheta)/theta^2*wBracket + (theta-sinTheta)/theta^3*wBracket*wBracket;
					v = V\p;
				end
				phibrac(1:3,4) = v;
			end
		end
		
		%%
		function [w,flag] = reparam(w)
			% Reparametrizes the exponential map to avoid singularities.
			% This is useful when using w to parametrize a rotation matrix.
			%
			% If |w| is close to pi, we replace w by (1 - 2*pi/|w|)*w,
			% which is an equivalent rotation, but with better derivatives.
			% https://www.cs.cmu.edu/~spiff/moedit99/expmap.pdf
			flag = false;
			wnorm = norm(w);
			while wnorm > 1.5*pi
				flag = true;
				a = (1-2*pi/wnorm);
				w = a*w;
				wnorm = abs(a*wnorm);
			end
		end
		
		%%
		function [ E ] = randE(  )
			% Creates a random transformation matrix
			[Q,R] = qr(randn(3)); %#ok<ASGLU>
			if det(Q) < 0
				% Negate the Z-axis
				Q(:,3) = -Q(:,3);
			end
			E = [Q, randn(3,1); 0 0 0 1];
		end
		
		%%
		function m = inertiaCuboid(whd,density)
			% Gets the diagonal inertia of a cuboid with (width, height, depth)
			if size(whd,1) < size(whd,2)
				whd = whd';
			end
			m = zeros(6,1);
			mass = density * prod(whd);
			m(1) = (1/12) * mass * whd([2,3])' * whd([2,3]);
			m(2) = (1/12) * mass * whd([3,1])' * whd([3,1]);
			m(3) = (1/12) * mass * whd([1,2])' * whd([1,2]);
			m(4) = mass;
			m(5) = mass;
			m(6) = mass;
		end
		
		%%
		function [xlines,ylines,zlines] = drawAxis(E,r)
			% Draws xyz axes at E
			if nargin < 2
				r = 1;
			end
			x = E(1:3,1);
			y = E(1:3,2);
			z = E(1:3,3);
			p = E(1:3,4);
			xlines = [p,p+r*x];
			ylines = [p,p+r*y];
			zlines = [p,p+r*z];
			if nargout == 0
				plot3(xlines(1,:),xlines(2,:),xlines(3,:),'r');
				plot3(ylines(1,:),ylines(2,:),ylines(3,:),'g');
				plot3(zlines(1,:),zlines(2,:),zlines(3,:),'b');
			end
		end
		
		%%
		function oneline = drawCuboid(E,whd,c)
			% Draws a cuboid at E with (width, height, depth)
			if nargin < 3
				c = 'k';
			end
			whd = whd/2;
			verts = ones(4,8);
			verts(1:3,1) = [-whd(1), -whd(2), -whd(3)]';
			verts(1:3,2) = [ whd(1), -whd(2), -whd(3)]';
			verts(1:3,3) = [ whd(1),  whd(2), -whd(3)]';
			verts(1:3,4) = [-whd(1),  whd(2), -whd(3)]';
			verts(1:3,5) = [-whd(1), -whd(2),  whd(3)]';
			verts(1:3,6) = [ whd(1), -whd(2),  whd(3)]';
			verts(1:3,7) = [ whd(1),  whd(2),  whd(3)]';
			verts(1:3,8) = [-whd(1),  whd(2),  whd(3)]';
			lines{1} = verts(:,[1,2,3,4,1]);
			lines{2} = verts(:,[5,6,7,8,5]);
			lines{3} = verts(:,[1,5]);
			lines{4} = verts(:,[2,6]);
			lines{5} = verts(:,[3,7]);
			lines{6} = verts(:,[4,8]);
			oneline = zeros(3,0);
			for j = 1 : length(lines)
				line = lines{j};
				for k = 1 : size(line,2)
					line(:,k) = E * line(:,k);
				end
				oneline = [oneline,nan(3,1),line(1:3,:)]; %#ok<AGROW>
			end
			if nargout == 0
				plot3(oneline(1,:),oneline(2,:),oneline(3,:),c);
			end
		end
		
		%%
		function [F,V] = patchCuboid(E,whd)
			% Gets the faces and vertices for a cuboid at E with (width, height, depth)
			whd = whd/2;
			verts = ones(4,8);
			verts(1:3,1) = [-whd(1), -whd(2), -whd(3)]';
			verts(1:3,2) = [ whd(1), -whd(2), -whd(3)]';
			verts(1:3,3) = [-whd(1),  whd(2), -whd(3)]';
			verts(1:3,4) = [ whd(1),  whd(2), -whd(3)]';
			verts(1:3,5) = [-whd(1), -whd(2),  whd(3)]';
			verts(1:3,6) = [ whd(1), -whd(2),  whd(3)]';
			verts(1:3,7) = [-whd(1),  whd(2),  whd(3)]';
			verts(1:3,8) = [ whd(1),  whd(2),  whd(3)]';
			verts = E*verts;
			V = verts(1:3,:)';
			F = [
				1 2 4 3
				8 7 5 6
				1 2 6 5
				8 7 3 4
				1 3 7 5
				8 6 2 4
				];
		end
	end
end
