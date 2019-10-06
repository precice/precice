import precice_future

class Interface:
    def __init__(self, solver_name, solver_process_index, solver_process_size):

        self.interface = precice_future.Interface(solver_name, solver_process_index, solver_process_size)

    def configure(self, configuration_file_name):
        self.interface.configure(configuration_file_name)

    def initialize(self):
        return self.interface.initialize()

    def initialize_data(self):
        self.interface.initialize_data()

    def advance(self, computed_timestep_length):
        return self.interface.advance(computed_timestep_length)

    def finalize(self):
        self.interface.finalize()

    def get_dimensions(self):
        return self.interface.get_dimensions()

    def is_coupling_ongoing(self):
        return self.interface.is_coupling_ongoing()

    def is_read_data_available(self):
        return self.interface.is_read_data_available()

    def is_write_data_required(self, computed_timestep_length):
        return self.interface.is_write_data_required(computed_timestep_length)

    def is_timestep_complete(self):
        return self.interface.is_timestep_complete()

    def is_action_required(self, action):
        return self.interface.is_action_required(action)

    def fulfilled_action(self, action):
        self.interface.fulfilled_action(action)

    def has_mesh(self, mesh_name):
        return self.interface.has_mesh(mesh_name)

    def get_mesh_id(self, mesh_name):
        return self.interface.get_mesh_id(mesh_name)

    def get_mesh_ids(self):
        return self.interface.get_mesh_ids()

    def has_data(self, data_name, mesh_id):
        return self.interface.has_data(data_name, mesh_id)

    def get_data_id(self, data_name, mesh_id):
        return self.interface.get_data_id(data_name, mesh_id)

    def set_mesh_vertices(self, mesh_id, size, positions, ids):
        dimensions = self.get_dimensions()
        assert(positions.size == size * dimensions)
        out_ids = self.interface.set_mesh_vertices(mesh_id, positions.reshape((size, dimensions)))
        ids[:] = out_ids[:]

    def set_mesh_vertex(self, mesh_id, position):
        return self.interface.set_mesh_vertex(mesh_id, position)

    def get_mesh_vertices(self, mesh_id, size, ids, positions):
        assert (positions.size == size * self.get_dimensions())
        out_positions = self.interface.get_mesh_vertices(mesh_id, ids)
        positions[:] = out_positions.flatten()[:]

    def get_mesh_vertex_size(self, mesh_id):
        return self.interface.get_mesh_vertex_size(mesh_id)

    def get_mesh_vertex_ids_from_positions(self, mesh_id, size, positions, ids):
        dimensions = self.get_dimensions()
        assert (positions.size == size * dimensions)
        out_ids = self.interface.get_mesh_vertex_ids_from_positions(mesh_id, positions.reshape((size, dimensions)))
        ids[:] = out_ids[:]

    def set_mesh_edge(self, mesh_id, first_vertex_id, second_vertex_id):
        return self.interface.set_mesh_edge(mesh_id, first_vertex_id, second_vertex_id)

    def set_mesh_triangle(self, mesh_id, first_edge_id, second_edge_id, third_edge_id):
        self.interface.set_mesh_triangle(mesh_id, first_edge_id, second_edge_id, third_edge_id)

    def set_mesh_triangle_with_edges(self, mesh_id, first_vertex_id, second_vertex_id, third_vertex_id):
        self.interface.set_mesh_triangle_with_edges(mesh_id, first_vertex_id, second_vertex_id, third_vertex_id)

    def set_mesh_quad(self, mesh_id, first_edge_id, second_edge_id, third_edge_id, fourth_edge_id):
        self.interface.set_mesh_quad(mesh_id, first_edge_id, second_edge_id, third_edge_id, fourth_edge_id)

    def set_mesh_quad_with_edges(self, mesh_id, first_vertex_id, second_vertex_id, third_vertex_id, fourth_vertex_id):
        self.interface.set_mesh_quad_with_edges(mesh_id, first_vertex_id, second_vertex_id, third_vertex_id, fourth_vertex_id)

    def map_read_data_to(self, to_mesh_id):
        self.interface.map_read_data_to(to_mesh_id)

    def map_write_data_from(self, from_mesh_id):
        self.interface.map_write_data_from(from_mesh_id)

    def write_block_vector_data(self, data_id, size, value_indices, values):
        dimensions = self.get_dimensions()
        assert (values.size == size * dimensions)
        assert (value_indices.size == size)
        self.interface.write_block_vector_data(data_id, value_indices, values.reshape((size, dimensions)))

    def write_vector_data(self, data_id, value_index, value):
        self.interface.write_vector_data(data_id, value_index, value)

    def write_block_scalar_data(self, data_id, size, value_indices, values):
        assert (values.size == size)
        assert (value_indices.size == size)
        self.interface.write_block_scalar_data(data_id, value_indices, values)

    def write_scalar_data(self, data_id, value_index, value):
        self.interface.write_scalar_data(data_id, value_index, value)

    def read_block_vector_data(self, data_id, size, value_indices, values):
        assert (values.size == size * self.get_dimensions())
        assert (value_indices.size == size)
        out_values = self.interface.read_block_vector_data(data_id, value_indices)
        values[:] = out_values.flatten()[:]

    def read_vector_data(self, data_id, value_index, value):
        out_value = self.interface.read_vector_data(data_id, value_index)
        value[:] = out_value[:]

    def read_block_scalar_data(self, data_id, size, value_indices, values):
        assert (values.size == size)
        assert (value_indices.size == size)
        out_values = self.interface.read_block_scalar_data(data_id, value_indices)
        values[:] = out_values[:]

    def read_scalar_data(self, data_id, value_index, value):
        out_value = self.interface.read_scalar_data(data_id, value_index)
        value[:] = out_value[:]
