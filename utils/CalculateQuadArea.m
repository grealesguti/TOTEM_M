function area = CalculateQuadArea(nodepos)
    % Numerically calculates the area of a quad element given the
    % nodal coordinates of its 4 vertices.

    % Define the vertices of the quad element
    vertices = nodepos(1:4, :);

    % Calculate the area of the quad using the shoelace formula
    x = vertices(:, 1);
    y = vertices(:, 2);

    area = 0.5 * abs(sum(x(1:end-1).*y(2:end) - x(2:end).*y(1:end-1)) + x(end)*y(1) - x(1)*y(end));
end