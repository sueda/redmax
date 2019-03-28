classdef JointPrismatic < redmax.Joint
	% https://github.com/junggon/gear/blob/master/src/gjoint_prismatic.cpp
	
	%%
	properties
		axis
		sides
	end
	
	%%
	methods
		%%
		function this = JointPrismatic(parent,body,axis)
			this = this@redmax.Joint(parent,body,1);
			this.axis = axis;
			this.sides = [0.5 0.5 0.5];
		end
		
		%%
		function setGeometry(this,sides)
			this.sides = sides;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function update_(this)
			this.Q = eye(4);
			this.Q(1:3,4) = this.axis*this.q(1);
			this.S(4:6,1) = this.axis;
		end
		
		%%
		function generateContacts_(this)
			% The contact points will be generated in the 'A' frame, which
			% is a frame defined with its Z-axis as the direction of
			% translation of this joint.
			S = eye(4);
			S(1:3,1:3) = diag(this.sides);
			E_ja = eye(4);
			z = [0 0 1]';
			angle = acos(max(-1.0, min(this.axis'*z, 1.0)));
			E_ja(1:3,1:3) = se3.aaToMat(cross(this.axis,z),angle);
			E_ia = this.body.E_ij * S * E_ja;
			x0 = -0.5;
			x1 =  0.5;
			y0 = -0.5;
			y1 =  0.5;
			for kz = 1 : 2
				z = (kz-1.5);
				tan_a = [0 0 1]';
				tan_i = E_ia(1:3,:)*[tan_a;0];
				tan_i = tan_i/norm(tan_i);
				for kx = 1 : 2
					for ky = 1 : 2
						if kx == 1 && ky == 1
							pos_a = [x0 y0 z]';
							nor1_a = [1 0 0]';
							nor2_a = [0 1 0]';
						elseif kx == 1 && ky == 2
							pos_a = [x0 y1 z]';
							nor1_a = [1 0 0]';
							nor2_a = [0 -1 0]';
						elseif kx == 2 && ky == 1
							pos_a = [x1 y0 z]';
							nor1_a = [-1 0 0]';
							nor2_a = [0 1 0]';
						else
							pos_a = [x1 y1 z]';
							nor1_a = [-1 0 0]';
							nor2_a = [0 -1 0]';
						end
						pos_i = E_ia(1:3,:)*[pos_a;1];
						nor1_i = E_ia(1:3,:)*[nor1_a;0];
						nor2_i = E_ia(1:3,:)*[nor2_a;0];
						nor1_i = nor1_i/norm(nor1_i);
						nor2_i = nor2_i/norm(nor2_i);
						this.contacts{end+1}.pos_i = pos_i;
						this.contacts{end}.nor_i = nor1_i;
						this.contacts{end}.tan_i = tan_i;
						this.contacts{end+1}.pos_i = pos_i;
						this.contacts{end}.nor_i = nor2_i;
						this.contacts{end}.tan_i = tan_i;
					end
				end
			end
			% Index for tangent matrix
			nt = redmax.Scene.countT();
			t = length(this.contacts);
			this.idxT = nt + (1:t);
			nt = nt + t;
			redmax.Scene.countT(nt);
		end
		
		%%
		function T = computeTangentMatrix_(this,T)
			for i = 1 : length(this.contacts)
				contact = this.contacts{i};
				tan_i = contact.tan_i;
				pos_i = contact.pos_i;
				T(this.idxT(i),this.body.idxM) = tan_i'*se3.Gamma(pos_i);
				if ~isempty(this.parent)
					E_pi = this.parent.body.E_iw*this.body.E_wi;
					tan_p = E_pi(1:3,:)*[tan_i;0];
					pos_p = E_pi(1:3,:)*[pos_i;1];
					T(this.idxT(i),this.parent.body.idxM) = -tan_p'*se3.Gamma(pos_p);
				end
			end
		end
		
		%%
		function [bl,bu,idx] = computeFrictionLimits_(this,mu,SPathresh,bl,bu,idx)
			% Combine normal and binormal, since they share the same
			% tangent.
			for i0 = 1 : 2 : length(this.contacts)
				i1 = i0 + 1;
				a0 = abs(this.contacts{i0}.a);
				a1 = abs(this.contacts{i1}.a);
				a = a0 + a1;
				bl(this.idxT(i0)) = -mu*a;
				bu(this.idxT(i0)) =  mu*a;
				if a > SPathresh
					idx(end+1) = this.idxT(i0); %#ok<AGROW>
				end
			end
		end
		
		%%
		function draw_(this)
			E_wj = this.body.E_wi*this.body.E_ij;
			[F,V] = se3.patchCuboid(E_wj,this.sides);
			patch('Faces',F,'Vertices',V,'FaceColor',this.body.color);
		end
	end
end
