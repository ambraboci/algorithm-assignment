import turtle
import math


class BigWeatherOptimizer:
    def __init__(self, num_dynos, possible_bonds, bucket_cost, bond_cost):

        self.num_dynos = num_dynos
        self.bucket_cost = bucket_cost
        self.bond_cost = bond_cost

        # Build the adjacency list representation of the graph
        self.graph = {}
        for i in range(1, num_dynos + 1):
            self.graph[i] = []

        for bond in possible_bonds:
            dyno1, dyno2 = bond
            self.graph[dyno1].append(dyno2)
            self.graph[dyno2].append(dyno1)

        # Store the possible bonds for reference
        self.possible_bonds = possible_bonds

        self.bucket_dynos = []
        self.selected_bonds = []

        self.num_optimal_solutions = 1

        self.optimal_configs = []

    def find_connected_components(self):

        visited = set()
        components = []

        for dyno in range(1, self.num_dynos + 1):
            if dyno not in visited:
                component = []
                queue = [dyno]
                visited.add(dyno)

                while queue:
                    current = queue.pop(0)
                    component.append(current)

                    for neighbor in self.graph[current]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)

                components.append(component)

        return components

    def find_mst(self, component):

        if not component:
            return []

        start = component[0]
        visited = {start}
        mst_edges = []

        pq = [(self.bond_cost, start, neighbor) for neighbor in self.graph[start]]

        while pq and len(visited) < len(component):

            min_idx = 0
            for i in range(1, len(pq)):
                if pq[i][0] < pq[min_idx][0]:
                    min_idx = i

            cost, frm, to = pq.pop(min_idx)

            if to not in visited:
                visited.add(to)
                edge = (min(frm, to), max(frm, to))
                mst_edges.append(edge)

                for neighbor in self.graph[to]:
                    if neighbor not in visited:
                        pq.append((self.bond_cost, to, neighbor))

        return mst_edges

    def minimum_dominating_set(self, component):

        uncovered = set(component)
        bucket_dynos = []
        bonds = []

        while uncovered:
            best_node = None
            best_coverage = -1

            for node in component:
                coverage = 1 if node in uncovered else 0
                for neighbor in self.graph[node]:
                    if neighbor in uncovered:
                        coverage += 1

                if coverage > best_coverage:
                    best_coverage = coverage
                    best_node = node

            if best_node is None:
                break

            bucket_dynos.append(best_node)

            uncovered.discard(best_node)
            for neighbor in self.graph[best_node]:
                uncovered.discard(neighbor)
                if neighbor in component and neighbor not in bucket_dynos:
                    bond = (min(best_node, neighbor), max(best_node, neighbor))
                    if bond not in bonds:
                        bonds.append(bond)

        # Calculate cost
        cost = len(bucket_dynos) * self.bucket_cost + len(bonds) * self.bond_cost

        return cost, bucket_dynos, bonds

    def find_all_optimal_configurations(self, component, target_cost):

        optimal_configs = []

        # 1. Consider all-buckets configuration
        all_buckets_cost = len(component) * self.bucket_cost
        if all_buckets_cost == target_cost:
            optimal_configs.append((component.copy(), []))

        # 2. Consider MST-based configurations with each node as a potential bucket host
        mst_bonds = self.find_mst(component)

        for central_bucket in component:
            centralized_cost = self.bucket_cost + len(mst_bonds) * self.bond_cost
            if centralized_cost == target_cost:
                optimal_configs.append(([central_bucket], mst_bonds))

        # 3. Try different combinations of bucket placements with varying numbers of buckets
        if len(component) <= 10:

            for num_buckets in range(1, len(component)):
                self._generate_bucket_combinations(component, num_buckets, [], 0, optimal_configs, target_cost)
        else:

            nodes_by_degree = []
            for node in component:
                degree = len(self.graph[node])
                nodes_by_degree.append((degree, node))

            nodes_by_degree.sort(reverse=True)

            for num_buckets in range(1, min(10, len(component))):
                top_nodes = [node for _, node in nodes_by_degree[:min(20, len(component))]]
                self._generate_bucket_combinations(top_nodes, num_buckets, [], 0, optimal_configs, target_cost)

        unique_configs = []
        for buckets, bonds in optimal_configs:
            bucket_set = set(buckets)
            bond_set = set(bonds)

            is_unique = True
            for unique_buckets, unique_bonds in unique_configs:
                if bucket_set == set(unique_buckets) and bond_set == set(unique_bonds):
                    is_unique = False
                    break

            if is_unique:
                unique_configs.append((buckets, bonds))

        return unique_configs

    def _generate_bucket_combinations(self, nodes, num_buckets, current_buckets, start_idx, optimal_configs,
                                      target_cost):

        if len(current_buckets) == num_buckets:
            bucket_set = set(current_buckets)
            bonds = []
            covered_nodes = set(current_buckets)

            for node in nodes:
                if node not in covered_nodes:
                    for bucket_node in current_buckets:
                        if node in self.graph[bucket_node]:
                            bond = (min(node, bucket_node), max(node, bucket_node))
                            if bond not in bonds:
                                bonds.append(bond)
                                covered_nodes.add(node)
                                break

            # If all nodes are covered, calculate the cost
            if len(covered_nodes) == len(nodes):
                cost = num_buckets * self.bucket_cost + len(bonds) * self.bond_cost
                if cost == target_cost:
                    optimal_configs.append((current_buckets.copy(), bonds))

            return

        for i in range(start_idx, len(nodes)):
            current_buckets.append(nodes[i])
            self._generate_bucket_combinations(nodes, num_buckets, current_buckets, i + 1, optimal_configs, target_cost)
            current_buckets.pop()

    def optimize_component(self, component):

        if not component:
            return 0, [], []

        # Strategy 1: Place buckets on all dynos (no bonds)
        all_buckets_cost = len(component) * self.bucket_cost

        min_cost = all_buckets_cost
        bucket_dynos = component.copy()
        selected_bonds = []

        # Strategy 2: Use MST-based approach (which can be optimal regardless of cost relationship)
        # Try each node as a potential bucket host for the MST-based approach
        mst_bonds = self.find_mst(component)

        for central_bucket in component:
            centralized_cost = self.bucket_cost + len(mst_bonds) * self.bond_cost

            if centralized_cost < min_cost:
                min_cost = centralized_cost
                bucket_dynos = [central_bucket]
                selected_bonds = mst_bonds

        # Strategy 3: Try dominating set approach
        dom_cost, dom_buckets, dom_bonds = self.minimum_dominating_set(component)

        if dom_cost < min_cost:
            min_cost = dom_cost
            bucket_dynos = dom_buckets
            selected_bonds = dom_bonds

        # Strategy 4: For case when bucket_cost < bond_cost, consider strategic bucket placements
        if self.bucket_cost < self.bond_cost:

            nodes_by_degree = []
            for node in component:
                degree = len(self.graph[node])
                nodes_by_degree.append((degree, node))

            nodes_by_degree.sort(reverse=True)

            for num_buckets in range(1, len(component)):
                test_buckets = [node for _, node in nodes_by_degree[:num_buckets]]

                test_bonds = []
                covered_nodes = set(test_buckets)

                for node in component:
                    if node not in covered_nodes:

                        best_dist = float('inf')
                        best_bucket = None

                        for bucket_node in test_buckets:
                            if node in self.graph[bucket_node]:
                                best_dist = 1
                                best_bucket = bucket_node
                                break

                        if best_bucket is not None:
                            bond = (min(node, best_bucket), max(node, best_bucket))
                            if bond not in test_bonds:
                                test_bonds.append(bond)
                                covered_nodes.add(node)

                if len(covered_nodes) == len(component):
                    test_cost = len(test_buckets) * self.bucket_cost + len(test_bonds) * self.bond_cost

                    if test_cost < min_cost:
                        min_cost = test_cost
                        bucket_dynos = test_buckets
                        selected_bonds = test_bonds

        all_configs = self.find_all_optimal_configurations(component, min_cost)

        num_configs = len(all_configs)
        if num_configs > 0:
            self.num_optimal_solutions *= num_configs
            self.optimal_configs.append(all_configs)

        return min_cost, bucket_dynos, selected_bonds

    def find_minimum_cost(self):

        components = self.find_connected_components()
        total_cost = 0
        all_bucket_dynos = []
        all_selected_bonds = []

        self.num_optimal_solutions = 1
        self.optimal_configs = []

        for component in components:
            cost, bucket_dynos, selected_bonds = self.optimize_component(component)
            total_cost += cost
            all_bucket_dynos.extend(bucket_dynos)
            all_selected_bonds.extend(selected_bonds)

        # Store the final solution
        self.bucket_dynos = all_bucket_dynos
        self.selected_bonds = all_selected_bonds

        return total_cost

    def text_visualize_solution(self):

        visualization = "BigWeather Solution Visualization:\n"
        visualization += "-" * 40 + "\n"

        visualization += "Bucket-hosting dynos: " + ", ".join(map(str, self.bucket_dynos)) + "\n"

        visualization += "Selected bonds: " + ", ".join([f"({a},{b})" for a, b in self.selected_bonds]) + "\n"

        bucket_cost = len(self.bucket_dynos) * self.bucket_cost
        bond_cost = len(self.selected_bonds) * self.bond_cost
        total_cost = bucket_cost + bond_cost

        visualization += f"Bucket cost: {len(self.bucket_dynos)} × {self.bucket_cost} = {bucket_cost}\n"
        visualization += f"Bond cost: {len(self.selected_bonds)} × {self.bond_cost} = {bond_cost}\n"
        visualization += f"Total cost: {total_cost}\n"

        visualization += f"Number of distinct optimal solutions: {self.num_optimal_solutions}\n"

        if self.optimal_configs and len(self.optimal_configs) > 0:
            visualization += "\nSample alternative optimal solutions:\n"

            for i, configs in enumerate(self.optimal_configs):
                if len(configs) > 1:
                    visualization += f"\nComponent {i + 1} alternatives:\n"

                    for j, (buckets, bonds) in enumerate(configs[:min(3, len(configs))]):
                        visualization += f"  Alternative {j + 1}:\n"
                        visualization += f"    Buckets: {', '.join(map(str, buckets))}\n"
                        visualization += f"    Bonds: {', '.join([f'({a},{b})' for a, b in bonds])}\n"

        visualization += "-" * 40 + "\n"

        return visualization


def visualize_solution(num_dynos, possible_bonds, bucket_dynos, selected_bonds, bucket_cost, bond_cost):
    screen = turtle.Screen()
    screen.title("BigWeather Solution Visualization")
    screen.setup(800, 600)
    screen.bgcolor("white")

    t = turtle.Turtle()
    t.speed(0)
    t.hideturtle()

    graph = {}
    for i in range(1, num_dynos + 1):
        graph[i] = []

    for bond in possible_bonds:
        dyno1, dyno2 = bond
        graph[dyno1].append(dyno2)
        graph[dyno2].append(dyno1)

    def find_connected_components():
        visited = set()
        components = []

        for dyno in range(1, num_dynos + 1):
            if dyno not in visited:
                component = []
                queue = [dyno]
                visited.add(dyno)

                while queue:
                    current = queue.pop(0)
                    component.append(current)

                    for neighbor in graph[current]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)

                components.append(component)

        return components

    components = find_connected_components()

    positions = {}

    num_components = len(components)
    big_radius = 0

    if num_components == 1:
        component_centers = [(0, 0)]
        big_radius = 0
    else:
        big_radius = 200
        component_centers = []
        for i in range(num_components):
            angle = 2 * math.pi * i / num_components
            x = big_radius * math.cos(angle)
            y = big_radius * math.sin(angle)
            component_centers.append((x, y))

    for idx, component in enumerate(components):
        center_x, center_y = component_centers[idx]
        component_size = len(component)

        radius = 30 + 15 * min(component_size, 10)

        for i, dyno in enumerate(component):
            angle = 2 * math.pi * i / component_size
            x = center_x + radius * math.cos(angle)
            y = center_y + radius * math.sin(angle)
            positions[dyno] = (x, y)

    t.pensize(1)
    t.pencolor("gray")

    for bond in possible_bonds:
        dyno1, dyno2 = bond

        if (dyno1, dyno2) in selected_bonds or (dyno2, dyno1) in selected_bonds:
            continue

        x1, y1 = positions[dyno1]
        x2, y2 = positions[dyno2]

        t.penup()
        t.goto(x1, y1)

        distance = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        dx = (x2 - x1) / distance
        dy = (y2 - y1) / distance

        node_radius = 20
        start_x = x1 + dx * node_radius
        start_y = y1 + dy * node_radius
        end_x = x2 - dx * node_radius
        end_y = y2 - dy * node_radius

        t.goto(start_x, start_y)

        dash_length = 5
        gap_length = 5
        current_distance = 0

        while current_distance < distance - 2 * node_radius:

            t.pendown()
            new_x = start_x + dx * min(dash_length, distance - 2 * node_radius - current_distance)
            new_y = start_y + dy * min(dash_length, distance - 2 * node_radius - current_distance)
            t.goto(new_x, new_y)
            current_distance += dash_length

            if current_distance >= distance - 2 * node_radius:
                break

            t.penup()
            new_x = start_x + dx * min(current_distance + gap_length, distance - 2 * node_radius)
            new_y = start_y + dy * min(current_distance + gap_length, distance - 2 * node_radius)
            t.goto(new_x, new_y)
            current_distance += gap_length

    t.pensize(2)
    t.pencolor("blue")

    for bond in selected_bonds:
        dyno1, dyno2 = bond
        x1, y1 = positions[dyno1]
        x2, y2 = positions[dyno2]

        distance = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        dx = (x2 - x1) / distance
        dy = (y2 - y1) / distance

        node_radius = 20
        start_x = x1 + dx * node_radius
        start_y = y1 + dy * node_radius
        end_x = x2 - dx * node_radius
        end_y = y2 - dy * node_radius

        t.penup()
        t.goto(start_x, start_y)
        t.pendown()
        t.goto(end_x, end_y)

    node_radius = 20

    for dyno in range(1, num_dynos + 1):
        x, y = positions[dyno]

        t.penup()
        t.goto(x, y - node_radius)

        if dyno in bucket_dynos:
            t.fillcolor("green")
        else:
            t.fillcolor("lightblue")

        t.pencolor("black")
        t.pensize(1)
        t.begin_fill()
        t.circle(node_radius)
        t.end_fill()

        t.penup()
        t.goto(x, y - 6)
        t.pencolor("black")
        t.write(str(dyno), align="center", font=("Arial", 12, "bold"))

    t.penup()
    t.goto(-350, 250)
    t.pencolor("black")
    t.write("Legend:", font=("Arial", 14, "bold"))

    t.penup()
    t.goto(-350, 220)
    t.fillcolor("lightblue")
    t.begin_fill()
    t.circle(10)
    t.end_fill()
    t.penup()
    t.goto(-330, 215)
    t.write("Regular Dyno", font=("Arial", 12))

    t.penup()
    t.goto(-350, 190)
    t.fillcolor("green")
    t.begin_fill()
    t.circle(10)
    t.end_fill()
    t.penup()
    t.goto(-330, 185)
    t.write("Bucket-Hosting Dyno", font=("Arial", 12))

    t.penup()
    t.goto(-350, 160)
    t.pencolor("gray")
    t.pendown()
    t.goto(-330, 160)
    t.penup()
    t.goto(-330, 155)
    t.pencolor("black")
    t.write("Possible Bond", font=("Arial", 12))

    t.penup()
    t.goto(-350, 130)
    t.pencolor("blue")
    t.pensize(2)
    t.pendown()
    t.goto(-330, 130)
    t.penup()
    t.goto(-330, 125)
    t.pencolor("black")
    t.pensize(1)
    t.write("Selected Bond", font=("Arial", 12))

    bucket_cost_total = len(bucket_dynos) * bucket_cost
    bond_cost_total = len(selected_bonds) * bond_cost
    total_cost = bucket_cost_total + bond_cost_total

    t.penup()
    t.goto(200, 250)
    t.pencolor("black")
    t.write("Cost Summary:", font=("Arial", 14, "bold"))

    t.penup()
    t.goto(200, 220)
    t.write(f"Buckets: {len(bucket_dynos)} × {bucket_cost} = {bucket_cost_total}", font=("Arial", 12))

    t.penup()
    t.goto(200, 190)
    t.write(f"Bonds: {len(selected_bonds)} × {bond_cost} = {bond_cost_total}", font=("Arial", 12))

    t.penup()
    t.goto(200, 160)
    t.write(f"Total Cost: {total_cost}", font=("Arial", 12, "bold"))

    t.penup()
    t.goto(200, 130)
    t.write(f"Number of dynos: {num_dynos}", font=("Arial", 12))

    t.penup()
    t.goto(200, 100)
    t.write(f"Number of components: {len(components)}", font=("Arial", 12))

    screen.mainloop()


def read_input(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        params = list(map(int, lines[0].strip().split()))
        num_dynos, num_bonds, bucket_cost, bond_cost = params

        possible_bonds = []
        for i in range(1, min(num_bonds + 1, len(lines))):
            bond = list(map(int, lines[i].strip().split()))
            possible_bonds.append(tuple(bond))

        return num_dynos, possible_bonds, bucket_cost, bond_cost


def main():
    file_path = input("Enter the input file path: ")

    try:
        num_dynos, possible_bonds, bucket_cost, bond_cost = read_input(file_path)
        optimizer = BigWeatherOptimizer(num_dynos, possible_bonds, bucket_cost, bond_cost)
        min_cost = optimizer.find_minimum_cost()

        print(f"The minimum cost for the BigWeather system is: {min_cost}")

        num_solutions = optimizer.num_optimal_solutions
        print(f"Number of distinct optimal solutions: {num_solutions}")

        print(optimizer.text_visualize_solution())

        show_visualization = input("Would you like to see the graphical visualization? (yes/no): ")
        if show_visualization.lower() == 'yes':
            visualize_solution(
                num_dynos,
                possible_bonds,
                optimizer.bucket_dynos,
                optimizer.selected_bonds,
                bucket_cost,
                bond_cost
            )

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
