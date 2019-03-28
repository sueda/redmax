classdef BodyMeshObj < redmax.Body
	%BodyMeshObj A body loaded from an .obj file
	
	%%
	properties
		filename % obj filename
		V        % vertices
		F        % faces
		E_oi     % Where the inertial frame is wrt the obj frame
		E_io     % Where the obj frame is wrt the inertial frame
	end
	
	%%
	methods
		%%
		function this = BodyMeshObj(density,filename)
			this = this@redmax.Body(density);
			this.filename = filename;
			this.E_io = eye(4);
			this.processMesh();
		end
		
		%%
		function processMesh(this)
			[this.V,this.F] = redmax.BodyMeshObj.readOBJ(this.filename);
			
			% The density argument is optional. We can easily multiply the resulting
			% inertia (element-by-element) by any density to get the same result.
			
			% trisurf(F',V(1,:),V(2,:),V(3,:),'FaceColor',[0.26,0.33,1.0 ]);
			% axis equal
			% xlabel('X'); ylabel('Y'); zlabel('Z');
			
			% Call VolInt
			[T0,T1,T2,TP] = redmax.BodyMeshObj.VolumeIntegration(this.V,this.F);
			
			density = this.density;
			mass = density * T0;
			
			%/* compute center of mass */
			r(1) = T1(1) / T0;
			r(2) = T1(2) / T0;
			r(3) = T1(3) / T0;
			
			%/* compute inertia tensor */
			J(1,1) = density * (T2(2) + T2(3));
			J(2,2) = density * (T2(3) + T2(1));
			J(3,3) = density * (T2(1) + T2(2));
			J(1,2) = - density * TP(1);
			J(2,3) = - density * TP(2);
			J(3,1) = - density * TP(3);
			J(2,1) = J(1,2);
			J(3,2) = J(2,3);
			J(1,3) = J(3,1);
			
			%/* translate inertia tensor to center of mass */
			J(1,1) = J(1,1) - mass * (r(2)*r(2) + r(3)*r(3));
			J(2,2) = J(2,2) - mass * (r(3)*r(3) + r(1)*r(1));
			J(3,3) = J(3,3) - mass * (r(1)*r(1) + r(2)*r(2));
			J(2,1) = J(2,1) + mass * r(1) * r(2);
			J(3,2) = J(3,2) + mass * r(2) * r(3);
			J(1,3) = J(1,3) + mass * r(3) * r(1);
			J(1,2) = J(2,1);
			J(2,3) = J(3,2);
			J(3,1) = J(1,3);
			
			% fprintf('center of mass:  (%+12.6f,%+12.6f,%+12.6f)\n\n', r(1), r(2), r(3));
			%
			% fprintf('inertia tensor with origin at c.o.m. :\n');
			% fprintf('%+15.6f  %+15.6f  %+15.6f\n', J(1,1), J(1,2), J(1,3));
			% fprintf('%+15.6f  %+15.6f  %+15.6f\n', J(2,1), J(2,2), J(2,3));
			% fprintf('%+15.6f  %+15.6f  %+15.6f\n\n', J(3,1), J(3,2), J(3,3));
			
			% Eigenvalue decomposition to get aligned frame
			E = eye(4);
			this.I_i = zeros(1,6);
			[JV,JD] = eig(J);
			this.I_i(1:3) = diag(JD);
			this.I_i(4:6) = mass;
			E(1:3,1:3) = JV;
			E(1:3,4) = r;
			
			% Check for right-handedness
			x = E(1:3,1);
			y = E(1:3,2);
			z = E(1:3,3);
			if cross(x,y)'*z < 0.0
				E(1:3,3) = -z;
			end
			this.E_oi = E;
			this.E_io = se3.inv(E);
			
			% alpha(0.3);
			% s = 4;
			% hold on;
			% plot3([r(1),r(1)+s*x(1)], [r(2),r(2)+s*x(2)], [r(3),r(3)+s*x(3)], 'r');
			% plot3([r(1),r(1)+s*y(1)], [r(2),r(2)+s*y(2)], [r(3),r(3)+s*y(3)], 'g');
			% plot3([r(1),r(1)+s*z(1)], [r(2),r(2)+s*z(2)], [r(3),r(3)+s*z(3)], 'b');
			% hold off;
			% axis equal;
		end
		
		%%
		function computeInertia_(this) %#ok<MANU>
			% Computes inertia at body
			% Should be already set!
		end
		
		%%
		function draw_(this)
			se3.drawAxis(this.E_wi,1);
			nverts = size(this.V,2);
			V_ = this.E_wi * this.E_io * [this.V;ones(1,nverts)];
			trisurf(this.F',V_(1,:),V_(2,:),V_(3,:),'FaceColor',this.color);
		end
		
		%%
		function s = getAxisSize(this) %#ok<MANU>
			s = 1;
		end
	end
	
	%%
	methods (Static)
		%%
		function [V,F] = readOBJ(filename)
			% http://www.alecjacobson.com/weblog/?p=917
			%
			% Reads a .obj mesh file and outputs the vertex and face list
			% assumes a 3D triangle mesh and ignores everything but:
			% v x y z and f i j k lines
			% Input:
			%  filename  string of obj file's path
			%
			% Output:
			%  V  3 x number of vertices array of vertex positions
			%  F  3 x number of faces array of face indices
			%
			V = zeros(3,0);
			F = zeros(3,0);
			vertex_index = 1;
			face_index = 1;
			fid = fopen(filename,'rt');
			line = fgets(fid);
			while ischar(line)
				vertex = sscanf(line,'v %f %f %f');
				f_v = sscanf(line,'f %d %d %d');
				f_vt = sscanf(line, 'f %d/%d %d/%d %d/%d');
				f_vn = sscanf(line, 'f %d//%d %d//%d %d//%d');
				f_vtn = sscanf(line,'f %d/%d/%d %d/%d/%d %d/%d/%d');
				
				% see if line is vertex command if so add to vertices
				if(length(vertex)>=3)
					V(:,vertex_index) = vertex;
					vertex_index = vertex_index+1;
					% see if line is a face (v) command if so add to faces
				elseif(length(f_v)>=3)
					F(:,face_index) = f_v;
					face_index = face_index+1;
					% see if line is a face (vt) command if so add to faces
				elseif(length(f_vt)>=6)
					F(:,face_index) = f_vt(1:2:end);
					face_index = face_index+1;
					% see if line is a face (vn) command if so add to faces
				elseif(length(f_vn)>=6)
					F(:,face_index) = f_vn(1:2:end);
					face_index = face_index+1;
					% see if line is a face (vtn) command if so add to faces
				elseif(length(f_vtn)>=9)
					% remove normal and texture indices
					f_vtn = f_vtn(1:3:end);
					F(:,face_index) = f_vtn;
					face_index = face_index+1;
				else
					%fprintf('Ignored: %s',line);
				end
				
				line = fgets(fid);
			end
			fclose(fid);
		end
		
		%%
		% [ T0 Tx Ty Tz Txx Tyy Tzz Txy Tyz Tzx ] = VolumeIntegration(NodePositions,Triangles,Normals)
		% this function calculates the volume from an assumed closed surface mesh defined by
		% Triangles and Positions (Nodepositions).  The Triangles list must have a counterclockwise numbering
		% to ensure the outward normal, and the Normals must be outward.  Normals are assumed to be of
		% unit magnitude.  As this is matlab, the first index is 1, not 0.
		%
		% 	/*******************************************************
		%       *                                                      *
		% 	*  volInt.m - Matlab Port of volInt.c by Brian Mirtich *
		% 	*                                                      *
		%       *  Ported to Matlab by Andrew Gosline                  *
		%       *  www.ece.ubc.ca/~andrewg/                            *
		%       *                                                      *
		% 	*  This code computes volume integrals needed for      *
		% 	*  determining mass properties of polyhedral bodies.   *
		% 	*                                                      *
		% 	*  For more information, see the paper:                *
		% 	*                                                      *
		% 	*  Brian Mirtich, "Fast and Accurate Computation of    *
		% 	*  Polyhedral Mass Properties," journal of graphics    *
		% 	*  tools, volume 1, number 1, 1996.                    *
		% 	*                                                      *
		% 	*  This source code is public domain, and may be used  *
		% 	*  in any way, shape or form, free of charge.          *
		% 	*                                                      *
		% 	*  Copyright 1995 by Brian Mirtich                     *
		% 	*                                                      *
		% 	*  mirtich@cs.berkeley.edu                             *
		% 	*  http://www.cs.berkeley.edu/~mirtich                 *
		%       *                                                      *
		% 	*******************************************************/
		%
		% /*
		% 	Revision history
		%
		% 	26 Jan 1996	Program creation.
		%
		% 	 3 Aug 1996	Corrected bug arising when polyhedron density
		% 			is not 1.0.  Changes confined to function main().
		% 			Thanks to Zoran Popovic for catching this one.
		%
		% 	27 May 1997     Corrected sign error in translation of inertia
		% 	                product terms to center of mass frame.  Changes
		% 			confined to function main().  Thanks to
		% 			Chris Hecker.
		% */
		%%%%%%
		% 2018
		%
		% Removed the Normals argument.
		%%%%%%
		function [ T0, T1, T2, TP ] = VolumeIntegration( V, F )
			Xn = V';
			Triangles = F';
			T0 = 0; Tx = 0; Ty = 0; Tz = 0; Txx = 0; Tyy = 0; Tzz = 0; Txy = 0; Tyz = 0; Tzx = 0; %#ok<NASGU>
			T1 = zeros(3,1); T2 = zeros(3,1); TP = zeros(3,1);
			for i = 1:length(Triangles)
				% Compute face normal
				tri = Triangles(i,:);
				v0 = Xn(tri(1),:);
				v1 = Xn(tri(2),:);
				v2 = Xn(tri(3),:);
				d10 = v1 - v0;
				d20 = v2 - v0;
				normal = cross(d10,d20);
				Normal = normal / norm(normal);
				if norm(normal) < 1e-9
					fprintf('skipping bad triangle %d\n', i);
					continue;
				end

				nx = abs(Normal(1)); ny = abs(Normal(2)); nz = abs(Normal(3));
				if (nx > ny) && (nx > nz)
					C = 0;
				elseif (ny > nz)
					C = 1;
				else
					C = 2;
				end
				
				%     Cycle the Indicies so always in numerical order
				A = mod(C+1,3);
				B = mod(A+1,3);
				A = A+1; B = B+1; C = C+1;
				
				%     calculate the offset
				w = -Normal(1)*Xn(Triangles(i,1),1) - Normal(2)*Xn(Triangles(i,1),2) - Normal(3)*Xn(Triangles(i,1),3);
				
				%     This part is essentially just the ComputeProjectionIntegrals part
				Pa = 0; Pb = 0; P1 = 0; Paa = 0; Pab = 0; Pbb = 0; Paaa = 0; Paab = 0; Pabb = 0; Pbbb = 0;
				
				for j = 1:3
					a0 = Xn(Triangles(i,j),A);
					b0 = Xn(Triangles(i,j),B);
					a1 = Xn(Triangles(i, mod(j,3)+1),A);
					b1 = Xn(Triangles(i, mod(j,3)+1),B);
					da = a1 - a0;
					db = b1 - b0;
					a0_2 = a0 * a0;
					a0_3 = a0_2 * a0;
					a0_4 = a0_3 * a0;
					b0_2 = b0 * b0;
					b0_3 = b0_2 * b0;
					b0_4 = b0_3 * b0;
					a1_2 = a1 * a1;
					a1_3 = a1_2 * a1;
					b1_2 = b1 * b1;
					b1_3 = b1_2 * b1;
					
					C1 = a1 + a0;
					Ca = a1*C1 + a0_2;
					Caa = a1 * Ca + a0_3;
					Caaa = a1*Caa + a0_4;
					Cb = b1 * (b1 + b0) + b0_2;
					Cbb = b1*Cb + b0_3;
					Cbbb = b1*Cbb + b0_4;
					Cab = 3*a1_2 + 2*a1*a0 + a0_2;
					Kab = a1_2 +2*a1*a0 + 3*a0_2;
					Caab = a0*Cab + 4*a1_3;
					Kaab = a1*Kab + 4*a0_3;
					Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
					Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;
					
					P1 = P1 + (db * C1);
					Pa = Pa + (db * Ca);
					Paa = Paa + db*Caa;
					Paaa = Paaa + db*Caaa;
					Pb = Pb + (da * Cb);
					Pbb = Pbb + da*Cbb;
					Pbbb = Pbbb + da*Cbbb;
					Pab = Pab + db*(b1*Cab + b0*Kab);
					Paab = Paab + db*(b1*Caab + b0*Kaab);
					Pabb = Pabb + da*(a1*Cabb + a0*Kabb);
				end
				
				P1 = P1/2.0;
				Pa = Pa/6.0;
				Paa = Paa/12.0;
				Paaa = Paaa/20.0;
				Pb = Pb/-6.0;
				Pbb = Pbb/-12.0;
				Pbbb = Pbbb/-20.0;
				Pab = Pab/24.0;
				Paab = Paab/60.0;
				Pabb = Pabb/-60.0;
				%disp(sprintf('%e %e %e %e %e %e %e %e %e %e \n',P1,Pa,Paa,Paaa,Pb,Pbb,Pbbb,Pab,Paab,Pabb));
				
				%     This Part is just the Compute Face Integrals part
				k1 = 1 / Normal(C);
				k2 = k1 * k1;
				k3 = k2 * k1;
				k4 = k3 * k1;
				
				Fa = k1 * Pa;
				Fb = k1 * Pb;
				Fc = -k2 * (Normal(A)*Pa + Normal(B)*Pb + w*P1);
				
				Faa = k1 * Paa;
				Fbb = k1 * Pbb;
				Fcc = k3 * ((Normal(A))^2*Paa + 2*Normal(A)*Normal(B)*Pab + (Normal(B))^2*Pbb + ...
					w*(2*(Normal(A)*Pa + Normal(B)*Pb) + w*P1));
				
				Faaa = k1 * Paaa;
				Fbbb = k1 * Pbbb;
				Fccc = -k4 * (Normal(A)^3*Paaa + 3*Normal(A)^2 * Normal(B)*Paab ...
					+ 3*Normal(A)*(Normal(B))^2*Pabb + (Normal(B))^3*Pbbb ...
					+ 3*w*((Normal(A))^2*Paa + 2*Normal(A)*Normal(B)*Pab + (Normal(B))^2*Pbb) ...
					+ w*w*(3*(Normal(A)*Pa + Normal(B)*Pb) + w*P1));
				
				Faab = k1 * Paab;
				Fbbc = -k2 * (Normal(A)*Pabb + Normal(B)*Pbbb + w*Pbb);
				Fcca = k3 * ((Normal(A))^2*Paaa + 2*Normal(A)*Normal(B)*Paab + (Normal(B))^2*Pabb ...
					+ w*(2*(Normal(A)*Paa + Normal(B)*Pab) + w*Pa));
				
				if A == 1
					Part = Fa;
				elseif B == 1
					Part = Fb;
				else
					Part = Fc;
				end
				
				T0 = T0 + Normal(1)*Part;
				T1(A) = T1(A) + Normal(A) * Faa;
				T1(B) = T1(B) + Normal(B) * Fbb;
				T1(C) = T1(C) + Normal(C) * Fcc;
				T2(A) = T2(A) + Normal(A) * Faaa;
				T2(B) = T2(B) + Normal(B) * Fbbb;
				T2(C) = T2(C) + Normal(C) * Fccc;
				TP(A) = TP(A) + Normal(A) * Faab;
				TP(B) = TP(B) + Normal(B) * Fbbc;
				TP(C) = TP(C) + Normal(C) * Fcca;
			end
			
			T1(1) = T1(1)/2; T1(2) = T1(2)/2; T1(3) = T1(3)/2;
			T2(1) = T2(1)/3; T2(2) = T2(2)/3; T2(3) = T2(3)/3;
			TP(1) = TP(1)/2; TP(2) = TP(2)/2; TP(3) = TP(3)/2;
		end
	end
end
