random_rotation(q::Node) = AngleAxis(q, rand(3)...) #make a random rotation function parameterized by variable q
random_rotation(q::Num) = random_rotation(Node(q))