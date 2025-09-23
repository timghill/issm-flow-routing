function edges = reconstruct_edges(md)
%Reconstruct edges
% {{{
% tic
%Maximum number of edges
maxnbf = 3*md.mesh.numberofelements;
%Initialize intermediaries
edges = zeros(maxnbf,3);
exchange = zeros(maxnbf,1);
%Chaining algorithm
head_minv = -1*ones(md.mesh.numberofvertices,1);
next_face = zeros(maxnbf,1);
%Initialize number of faces
nbf = 0;
for i=1:md.mesh.numberofelements
	for j=1:3
		%Get the two indices of the edge number j of the ith triangle
		v1 = md.mesh.elements(i,j);
		if(j==3)
			v2 = md.mesh.elements(i,1);
		else
			v2 = md.mesh.elements(i,j+1);
		end
		%sort
		if(v2<v1)
			v3=v2; v2=v1; v1=v3;
		end
		%This edge a priori has not been processed yet
		exists = false;
		%Go through all processed faces connected to v1 and check whether we have seen this edge yet
		e=head_minv(v1);
		while(e~=-1)
			if(edges(e,2)==v2)
				exists = true;
				break;
			end
			e=next_face(e);
		end
		%If this edge is new, add it to the lists
		if(~exists)
			%Update edges
			edges(nbf+1,1) = v1; %vertex 1
			edges(nbf+1,2) = v2; %vertex 2
			edges(nbf+1,3) = i;  %element 1 (ignore second one)
			if(v1~=md.mesh.elements(i,j)) exchange(nbf+1)=1; end
			%Update chain
			next_face(nbf+1) = head_minv(v1);
			head_minv(v1)    = nbf+1;
			%Increase number of faces
			nbf=nbf+1;
		end
	end
end

edges = edges(1:nbf,:);
pos = find(exchange);
v3 = edges(pos,1);
edges(pos,1) = edges(pos,2);
edges(pos,2) = v3;
% toc
% }}}

%Change edges formatting so that plot looks ok
myedges = [edges(:,1:2) edges(:,1)]';