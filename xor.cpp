#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <memory>
#include <map>

#include <iostream>

std::random_device rd;
std::mt19937 gen(rd()); 

// nerual network config
double c[3] = { 0.5, 0.5, 0.3 };
double weight_delta = 0.01;

// family config
int survivors = 1;
int elitism = 3;
double family_threshold = 3.0;

// mutation config
double mutate_chance = 0.7;
double change_weight_chance = 0.20;
double adjust_weight_chance = 0.65;
double add_connection_chance = 0.07;
double split_connection_chance = 0.03;
double toggle_connection_chance = 0.05;


template <typename Iter>
auto& random_element(Iter begin, Iter end)
{
	int d = static_cast<int>(std::distance(begin, end)) - 1;
	int a = std::uniform_int_distribution<int>{ 0, d }(gen);
	std::advance(begin, a);
	return *begin;
}

auto& random_element(auto&& c)
{
	return random_element(c.begin(), c.end());
}

double random_number(double min = -1.0, double max = 1.0)
{
	return std::uniform_real_distribution<double>{ min, std::nextafter(max, std::numeric_limits<double>::max()) }(gen);
}

double random_probability()
{
	return std::uniform_real_distribution<double>{ 0.0, 1.0 }(gen);
}


long combine(long x, long y)
{
	return x << 32 | y;
}

class Node
{
	static int recent_innovation;
	static std::unordered_map<long, int> discovered_nodes;
public:
	Node() : id{ -1 }
	{}
	int id;
	double x;
	double y;

	std::unordered_set<int> connections;

	Node(double x, double y, int id) : id{ id }, x{ x }, y{ y }
	{
		if (recent_innovation <= id)
			recent_innovation = id + 1;
	}

	Node(Node const& in, Node const& out) : x{ (in.x + out.x) / 2 }, y{ (in.y + out.y) / 2 * random_number(0.9, 1.1) }
	{
		auto comb = combine(in.id, out.id);
		if (!discovered_nodes.contains(comb))
			discovered_nodes[comb] = recent_innovation++;
		id = discovered_nodes[comb];
	}

	void add_connection(int node_id)
	{
		connections.insert(node_id);
	}

	bool operator<(Node const& n) const
	{
		return x < n.x;
	}

	bool connected_to(Node const& n) const
	{
		return connections.contains(n.id);
	}
};

int Node::recent_innovation = 0;
std::unordered_map<long, int> Node::discovered_nodes;

class Connection
{
	static int recent_innovation;
	static std::unordered_map<long, int> discovered_connections;
public:
	Connection() : id{ -1 }
	{}

	int id;

	double weight;
	int in;
	int out;
	bool enabled;

	Connection(Node& in, Node& out, double weight = 0.0, bool enabled = true) : in{ in.id }, out{ out.id }, weight{ weight }, enabled{ enabled }
	{
		auto comb = combine(in.id, out.id);
		if (!discovered_connections.contains(comb))
			discovered_connections[comb] = recent_innovation++;
		id = discovered_connections[comb]; 
		out.add_connection(in.id);
	}

	static int id_for(int in, int out)
	{
		auto com = combine(in, out);
		return discovered_connections.contains(com) ? discovered_connections[com] : -1; 
	}
};

int Connection::recent_innovation = 0;

std::unordered_map<long, int> Connection::discovered_connections;

template <typename T>
class Id_Map
{
	std::map<int, T> data;
	using Base_Iterator = typename std::map<int, T>::iterator;
	using Base_Const_Iterator = typename std::map<int, T>::const_iterator;
public:
	class iterator 
	{
		Base_Iterator iter;
	public:
		using difference_type = Base_Iterator::difference_type;
		using value_type = Base_Iterator::value_type;
		using pointer = Base_Iterator::pointer;
		using reference =  Base_Iterator::reference;
		using iterator_category =  Base_Iterator::iterator_category;

		iterator(Base_Iterator iter) : iter{ iter }
		{}
		auto& operator*()
		{
			return iter->second;
		}
		auto* operator->()
		{
			return &iter->second;
		}
		void operator++()
		{
			++iter;
		}
		void operator--()
		{
			--iter;
		}
		iterator operator+(int amount)
		{
			auto i = iter;
			std::advance(i, amount);
			return i;
		}
		bool operator!=(iterator other) const
		{
			return iter != other.iter;
		}
	};
	class const_iterator
	{
		Base_Const_Iterator iter;
	public:
		using difference_type = Base_Const_Iterator::difference_type;

		const_iterator(Base_Const_Iterator iter) : iter{ iter }
		{}
		auto const& operator*()
		{
			return iter->second;
		}
		auto const* operator->()
		{
			return &iter->second;
		}
		void operator++()
		{
			++iter;
		}
		bool operator!=(const_iterator other) const
		{
			return iter != other.iter;
		}
	};
	iterator begin() { return data.begin(); }
	const_iterator begin() const { return data.begin(); }

	iterator end() { return data.end(); }
	const_iterator end() const { return data.end(); }

	void insert(auto&& t) {	data.insert({t.id, t}); }
	void remove(auto&& t) { data.erase(t.id); }
	
	auto size() const {	return data.size(); }

	auto& operator[](int id) { return data[id]; }
	auto const& operator[](int id) const { return data.at(id); }

	bool contains(T const& t) const { return data.contains(t.id); }
};


//auto cmp_id = [](auto&& x, auto&& y){ return x.id < y.id; };

namespace Graphic
{
	template <auto, typename>
	class Network;
}

template <typename Type>
class Network
{
	Id_Map<Node> nodes;
	Id_Map<Connection> genes;

	std::vector<int> outputs;
	std::vector<int> inputs;

	double activation(double x) const
	{
		return 1.0 / (1.0 + std::exp(x));
	}

	template <typename T>
	using Gene_Types = std::array<T, 3>;

	enum class Gene_Type
	{
		matching = 0,
		disjoint = 1,
		excess = 2,
	};

	void detail_gene_types(Type const& other, auto&& found) const
	{
		auto my_iter = genes.begin();
		auto other_iter = other.genes.begin();

		auto my_end = genes.end();
		auto other_end = other.genes.end();

		while (my_iter != my_end && other_iter != other_end)
			if (my_iter->id == other_iter->id)
			{
				found(Gene_Type::matching, 0, *my_iter);
				found(Gene_Type::matching, 1, *other_iter);
				++my_iter;
				++other_iter;
			}
			else if (my_iter->id < other_iter->id)
			{
				found(Gene_Type::disjoint, 0, *my_iter);
				++my_iter;
			}
			else 
			{
				found(Gene_Type::disjoint, 1, *other_iter);
				++other_iter;
			}

		for (;my_iter != my_end; ++my_iter)
			found(Gene_Type::excess, 0, *my_iter);
		for (;other_iter != other_end; ++other_iter)
			found(Gene_Type::excess, 1, *other_iter);
	}

	std::tuple<int, int, double> distance_gene_types(Type const& other) const
	{
		std::tuple<double, int, int> r;
		int m = 0;
		double w = 0;
		detail_gene_types(other, [&](Gene_Type type, int parent, Connection const& c)
		{
			switch (type)
			{
			case Gene_Type::matching:
				m++;
				if (m & 1)
					std::get<0>(r) += std::abs(w - c.weight);
				else
					w = c.weight;
				return;
			case Gene_Type::disjoint:
				std::get<1>(r)++;
				return;
			case Gene_Type::excess:
				std::get<2>(r)++;
				return;
			}
		});
		std::get<0>(r) /= m;
		return r;
	}

	std::array<Gene_Types<std::vector<Connection>>, 2> gene_types(Type const& other) const
	{
		std::array<Gene_Types<std::vector<Connection>>, 2> r;
		detail_gene_types(other, [&](Gene_Type type, int parent, Connection const& c)
		{
			r[parent][static_cast<int>(type)].push_back(c);
		});
		return r;
	}
	double operator()(int out, std::unordered_map<int, double>& dic) const 
	{
		if (!dic.contains(out))
		{
			double r = 0;
			for (auto in : nodes[out].connections)
			{
				long id = Connection::id_for(in, out);
				if (genes[id].enabled)
					r += genes[id].weight * (*this)(in, dic);
			}
			dic[out] = activation(r);
		}
		return dic[out];
	}

	void derive_nodes()
	{
		std::unordered_set<int> existing_nodes;
		for (auto id : inputs)
			existing_nodes.insert(id);
		for (auto id : outputs)
			existing_nodes.insert(id);
		std::unordered_map<int, int> nodes_to_add;
		for (auto& gene : genes)
		{
			if (!existing_nodes.contains(gene.out))
				nodes_to_add[gene.out] = gene.in;
			if (!existing_nodes.contains(gene.in))
			{
				auto n = Node{ nodes[nodes_to_add[gene.in]], nodes[gene.out] };
				nodes_to_add.erase(gene.in);
				existing_nodes.insert(n.id);
				nodes.insert(n);
			}
		}
		for (auto& gene : genes)
			nodes[gene.out].add_connection(gene.in);
	}
public:
	void add_input(Node const& node)
	{
		inputs.push_back(node.id);
		nodes.insert(node);
	}

	void add_output(Node const& node)
	{
		outputs.push_back(node.id);
		nodes.insert(node);
	}

	std::vector<double> operator()(std::vector<double> const& input) const
	{
		std::unordered_map<int, double> d;
		for (int i = 0; i < inputs.size(); i++)
			d[inputs[i]] = input[i];
		std::vector<double> r(outputs.size());
		int c = 0;
		for (auto o : outputs)
			r[c++] = (*this)(o, d);
		return r;
	}

	Type crossover_with(Type const& other) const
	{
		auto kid = Type{};
		auto classified_genes = gene_types(other);
		auto&& [m, d, e] = classified_genes[0];
		auto&& n = classified_genes[1][0];

		for (int i = 0; i < m.size(); i++)
			kid.genes.insert(random_probability() < 0.5 ? m[i] : n[i]);
		for (auto& disjoint : d)
			kid.genes.insert(disjoint);
		for (auto& excess : e)
			kid.genes.insert(excess);

		kid.derive_nodes();

		return kid;
	}

	double distance_from(Type const& other) const
	{
		auto [w, d, e] = distance_gene_types(other);
		auto max = std::max(genes.size(), other.genes.size());
		return c[0] * d / max + c[1] * e / max + c[2] * w;
	}

	void randomly_split_connection()
	{
		auto& g = random_element(genes);
		auto n = Node(nodes[g.in], nodes[g.out]);
		if (!nodes.contains(n))
		{
			auto g1 = Connection(nodes[g.in], n, 1, true);
			auto g2 = Connection(n, nodes[g.out], g.weight, g.enabled);
			g.enabled = false;
			nodes.insert(n);
			genes.insert(g1);
			genes.insert(g2);
		}
	}

	void connect(int in, int out)
	{
		genes.insert(Connection{ nodes[in], nodes[out], random_number(), true });
	}

	void randomly_add_connection()
	{
		auto& n1 = random_element(nodes.begin(), nodes.end());
		auto& n2 = random_element(nodes.begin() + inputs.size(), nodes.end());
		if (n1 < n2)
			if (!n2.connected_to(n1))
				genes.insert(Connection{ n1, n2, random_number(), true });
		else if (n2 < n1)
			if (!n1.connected_to(n2))
				genes.insert(Connection{ n2, n1, random_number(), true });
	}

	void randomly_change_weight()
	{
		auto& g = random_element(genes);
		g.weight = random_number();
	}

	void randomly_adjust_weight()
	{
		auto& g = random_element(genes);
		g.weight += random_number(-weight_delta, weight_delta);
		g.weight = std::clamp(g.weight, -1.0, 1.0);
	}

	void randomly_toggle_connection()
	{
		auto& g = random_element(genes);
		g.enabled = !g.enabled;
	}

	void print_genes()
	{
		std::cout << "\033[2J\033[H" << std::endl;
		for (auto& gene : genes)
			std::cout << gene.id << " " << gene.in << " " << gene.out << " " << gene.weight << " " << gene.enabled << std::endl;
	}

	void print_nodes()
	{
		std::cout << "\033[2J\033[H" << std::endl;
		for (auto& node : nodes)
			std::cout << node.id << " (" << node.x << ", " << node.y << ")" << std::endl;
	}

	template <auto action, typename T>
	friend class Graphic::Network;
};

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>

namespace Graphic
{
	namespace SDL 
	{
		static class SDL_LIBS
		{
		public:
			SDL_LIBS()
			{
				SDL_Init(SDL_INIT_EVERYTHING);
				TTF_Init();
			}
			~SDL_LIBS()
			{
				TTF_Quit();
				SDL_Quit();
			}
		} sdl_libs{};

		class Window
		{
			SDL_Window* window;
		public:
			Window(std::string const& title, int x = 0, int y = 0, int w = 800, int h = 600, int flags = 0) : window{ SDL_CreateWindow(title.c_str(), x, y, w, h, flags | SDL_WINDOW_RESIZABLE ) }
			{}

			~Window()
			{
				if (window)
					SDL_DestroyWindow(window);
			}

			friend class Renderer;
		};

		class Renderer
		{
			SDL_Renderer* renderer;
			TTF_Font* font;
			class Draw_Color
			{
				SDL_Renderer* renderer;
				SDL_Color org;
			public:
				Draw_Color(SDL_Renderer* renderer, SDL_Color color) : renderer{ renderer }
				{
					SDL_GetRenderDrawColor(renderer, &org.r, &org.g, &org.b, &org.a);
					SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
				}
				~Draw_Color()
				{
					SDL_SetRenderDrawColor(renderer, org.r, org.g, org.b, org.a);
				}
			};
			[[nodiscard]]
			Draw_Color change_draw_color(SDL_Color color)
			{
				return Draw_Color{ renderer, color };
			}
		public:
			Renderer(Window& window, std::string const& font_file) : renderer{ SDL_CreateRenderer(window.window, -1, 0) }, font{ TTF_OpenFont(font_file.c_str(), 18) }
			{}	
			~Renderer()
			{
				if (renderer)
					SDL_DestroyRenderer(renderer);
				if (font)
					TTF_CloseFont(font);
			}
			void draw_rectangle(int x, int y, int w, int h, SDL_Color color)
			{
				auto state = change_draw_color(color);
				SDL_Rect r{ x, y, w, h };
				SDL_RenderDrawRect(renderer, &r);
			}
			void draw_circle(int x, int y, int r, SDL_Color color)
			{
				auto state = change_draw_color(color);
				SDL_Point points[64];
				for (int i = 0; i < 64; i++)
				{
					points[i].x = x + r / 2 + r * std::cos(2.0 * 3.14 * i / 63.0);
					points[i].y = y + r / 2 + r * std::sin(2.0 * 3.14 * i / 63.0);
				}
				SDL_RenderDrawLines(renderer, &points[0], 64);
			}
			std::array<int, 2> size_text(std::string const& text)
			{
				std::array<int, 2> r;
				TTF_SizeText(font, text.c_str(), &r[0], &r[1]);
				return r;
			}
			void draw_text(std::string const& text, int x, int y, SDL_Color color)
			{
				auto s = TTF_RenderText_Blended(font, text.c_str(), color);
				auto t = SDL_CreateTextureFromSurface(renderer, s);
				auto [w, h] = size_text(text);
				SDL_Rect r{ x, y, w, h };
				SDL_RenderCopy(renderer, t, NULL, &r);
				SDL_DestroyTexture(t);
				SDL_FreeSurface(s);
			}
			void draw_line(int x1, int y1, int x2, int y2, SDL_Color color)
			{
				auto state = change_draw_color(color);
				SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
			}
			void clear()
			{
				SDL_RenderClear(renderer);
			}
			void present()
			{
				SDL_RenderPresent(renderer);
			}
		};
	}

	class Drawable
	{
	public:
		virtual void draw(SDL::Renderer& r) = 0;
		virtual void on_click(int x, int y){}
	};

	class Window
	{
		std::vector<std::unique_ptr<Drawable>> controls;
		SDL::Window window;
		SDL::Renderer renderer;
		bool is_running;
		void handle_event(SDL_Event& e)
		{
			switch(e.type)
			{
			case SDL_QUIT:
				quit();
				break;
			case SDL_MOUSEBUTTONDOWN:
				for (auto& c : controls)
					c->on_click(e.button.x, e.button.y);
			}
		}
	public:
		Window(std::string const& title, std::string const& font_file) : window{ title }, renderer{ window, font_file }, is_running{ false }
		{}

		void add_control(auto&& c)
		{
			controls.emplace_back(std::make_unique<std::remove_reference_t<decltype(c)>>(std::forward<decltype(c)>(c)));
		}
		void draw()
		{
			renderer.clear();
			for (auto& control : controls)
				control->draw(renderer);
			renderer.present();
		}
		int run()
		{
			is_running = true;
			SDL_Event e;
			while (is_running)
			{
				while(SDL_PollEvent(&e))
					handle_event(e);
				draw();
				SDL_Delay(1000 / 30);
			}
			return 0;
		}
		void quit()
		{
			is_running = false;
		}
	};

	template <auto action, typename ... T>
	class Button : public Drawable
	{
		std::string text;
		std::tuple<T&...> n;
		int x, y, w, h;
		SDL_Color border, text_color;
	public:
		Button(std::string const& s, int x, int y, int w, int h, SDL_Color border, SDL_Color text_color, T& ... t) : 
			text{ s }, x{ x }, y{ y }, w{ w }, h{ h }, border{ border }, text_color{ text_color }, n{ std::tuple<T&...>{ t... } }
		{}
		void draw(SDL::Renderer& r) final
		{
			r.draw_rectangle(x, y, w, h, border);
			auto [tx, ty] = r.size_text(text);
			r.draw_text(text, x + (w - tx) / 2, y + (h - ty) / 2, text_color);
		}
		void on_click(int x, int y)
		{
			if (this->x < x && x < this->x + w && this->y < y && y < this->y + h)
				std::apply(action, n);
		}
	};

	template <auto action, typename T>
	class Network : public Drawable
	{
		T& t;
		int x, y, w, h;
		SDL_Color background, node, enabled_connection, disabled_connection, text;
	public:
		Network(int x, int y, int w, int h, SDL_Color background, SDL_Color node, SDL_Color enabled_connection, SDL_Color disabled_connection, SDL_Color text, T& n) :
			t{ n }, x{ x }, y{ y }, w{ w }, h{ h }, background{ background }, node{ node }, enabled_connection{ enabled_connection }, disabled_connection{ disabled_connection }, text{ text }
		{}

		void draw(SDL::Renderer& r) final
		{
			auto& n = action(t);
			r.draw_rectangle(x, y, w, h, background);
			for (auto const& gene : n.genes)
			{
				auto& n1 = n.nodes[gene.in];
				auto& n2 = n.nodes[gene.out];
				r.draw_line(x + n1.x + 15, y + n1.y + 15, x + n2.x + 15, y + n2.y + 15, gene.enabled ? enabled_connection : disabled_connection);
				auto label = std::to_string(gene.weight);
				auto [tx, ty] = r.size_text(label);
				r.draw_text(label, x + (n2.x + n1.x - tx) / 2, y + (n2.y + n1.y - ty) / 2, text);
			}
			for (auto const& node : n.nodes)
			{
				r.draw_circle(x + node.x, y + node.y, 30, this->node);
				auto label = std::to_string(node.id);
				auto [tx, ty] = r.size_text(label);
				r.draw_text(label, x + node.x + (30 - tx) / 2, y + node.y + (30 - ty) / 2, text);
			}
		}
	};
}

class Xor_Brain : public Network<Xor_Brain>
{
public:
	Xor_Brain()
	{
		add_input(Node{ 100, 100, 0 });
		add_input(Node{ 100, 200, 1 });
		add_input(Node{ 100, 300, 2 });
		add_output(Node{ 1500, 200, 3 });
	}
};

class Xor_Creature
{
	std::tuple<Xor_Creature const&, Xor_Creature const&> compare_to(Xor_Creature const& other) const
	{
		if (score >= other.score)
			return { *this, other };
		else
			return { other, *this };
	}
public:
	Xor_Brain brain;
	double score;
	Xor_Creature() : brain{}, score{ 0 }
	{
		brain.connect(0, 3);
		brain.connect(1, 3);
		brain.connect(2, 3);
	}

	Xor_Creature(Xor_Brain const& brain) : brain{ brain }, score{ 0 }
	{}

	Xor_Creature breed_with(Xor_Creature const& other) const
	{
		auto&& [b, w] = compare_to(other);
		return Xor_Creature{ b.brain.crossover_with(w.brain) };
	}

	bool is_related_to(Xor_Creature const& other) const
	{
		auto&& [b, w] = compare_to(other);
		return b.brain.distance_from(w.brain) < family_threshold;
	}

	void mutate()
	{
		auto p = random_probability(), cum = 0.0;
		if (p < (cum += change_weight_chance))
			brain.randomly_change_weight();
		else if (p < (cum += adjust_weight_chance))
			brain.randomly_adjust_weight();
		else if (p < (cum += toggle_connection_chance))
			brain.randomly_toggle_connection();
		else if (p < (cum += add_connection_chance))
			brain.randomly_add_connection();
		else if (p < (cum += split_connection_chance))
			brain.randomly_split_connection();
	}

	std::array<std::vector<double>, 4> results() const
	{
		return {
			brain({ 0, 0, 1 }),
			brain({ 0, 1, 1 }),
			brain({ 1, 0, 1 }),
			brain({ 1, 1, 1 })
		};
	}
	
	double evaluate() 
	{
		double correct[] = { 0, 1, 1, 0 };
		auto r = results();
		double d = 0;
		for (int i = 0; i < 4; i++)
			d += std::abs(correct[i] - r[i][0]);
		d = 4 - d;
		return score = d * d;
	}

	void print_results()
	{
		char const* pre[] = { "0,0", "0,1", "1,0", "1,1" };
		auto results = this->results();
		for (int i = 0; i < 4; i++)
			std::cout << pre[i] << ": " << results[i][0] << std::endl;
	}

	void print_score()
	{
		std::cout << score << std::endl;
	}
};


template <typename T>
std::tuple<T&, T&> random_weighted_pair(std::vector<T>& c, std::discrete_distribution<>& dis)
{
	if (c.size() == 1)
		return std::tuple<T&, T&>{ c.front(), c.front() };
	int p = dis(gen), q = dis(gen);
	if (p == q)
		if (q != 0)
			q--;
		else
			q++;
	return std::tuple<T&, T&>{ c[p], c[q] };
}

class Xor_Family
{
	std::vector<Xor_Creature> creatures;
public:
	Xor_Family() = default;
	Xor_Family(int size) : creatures(size)
	{}

	double score;
	
	auto begin() { return creatures.begin(); }
	auto begin() const { return creatures.begin(); }

	auto end() { return creatures.end(); }
	auto end() const { return creatures.end(); }

	void mutate()
	{
		for (auto& c : creatures)
			if (random_probability() < mutate_chance)
				c.mutate();
	}
	void breed(int size)
	{
		std::vector<Xor_Creature> next;

		int survivor_amount = survivors > size ? size : survivors;
		std::copy(creatures.begin(), creatures.begin() + survivor_amount, std::back_inserter(next));
		size -= survivor_amount;

		
		int pop_to_kill = creatures.size() / elitism;
		creatures.erase(creatures.end() - pop_to_kill, creatures.end()); 

		std::vector<double> weights;
		for (auto& c : creatures)
			weights.push_back(c.score);

		std::discrete_distribution<> dis{ weights.begin(), weights.end() };

		for (int i = 0; i < size; i++)
		{	
			auto&& [p, q] = random_weighted_pair(creatures, dis);	
			next.push_back(p.breed_with(q));
		}
			
		creatures = std::move(next);
	}
	double evaluate()
	{
		score = 0;
		for (auto& c : creatures)
			score += c.evaluate();
		std::sort(creatures.begin(), creatures.end(), [](auto const& x, auto const& y){ return x.score > y.score; });
		return score /= creatures.size();
	}
	bool is_related_to(Xor_Creature const& c) const
	{
		return creatures.size() ? head().is_related_to(c) : true;
	}
	bool is_related_to(Xor_Family const& f) const
	{
		return creatures.size() ? head().is_related_to(f.head()) : true;
	}

	Xor_Creature const& head() const
	{
		return creatures.front();
	}

	Xor_Creature& head()
	{
		return creatures.front();
	}

	void add(Xor_Family const& f)
	{
		for (auto& c : f)
			creatures.push_back(c);
	}

	void add(Xor_Creature const& c)
	{
		creatures.push_back(c);
	}

	std::vector<Xor_Creature> expell_misfits()
	{
		std::vector<Xor_Creature> misfits;
		creatures.erase(std::remove_if(creatures.begin(), creatures.end(), [&](auto& c)
		{ 
			bool r = !is_related_to(c);
			if (r)
				misfits.push_back(c);
			return r;
		}), creatures.end());
		return misfits;
	}

	void print_results()
	{
		for (auto& c : creatures)
			c.print_results();
	}

	void print_scores()
	{
		std::cout << "family score: " << score << std::endl;
		for (auto& c : creatures)
			c.print_score();
	}

	Xor_Creature& operator[](int index)
	{
		return creatures[index];
	}

	int size() const { return creatures.size(); }
};

class Xor_Specie
{
	std::vector<Xor_Family> families;
	int size;
	double score;
	Xor_Specie() = default;
public:
	Xor_Specie(int size) : families{ { Xor_Family{ size } } }, size{ size }
	{
	}

	auto begin() { return families.begin(); }
	auto begin() const { return families.begin(); }

	auto end() { return families.end(); }
	auto end() const { return families.end(); }
	
	void evolve(int generations)
	{
		evaluate();
		for (int i = 1; i < generations; i++)
		{
			breed();
			mutate();
			split();
			evaluate();
		}
	}

	void mutate()
	{
		for (auto& f : families)
			f.mutate();
	}

	void breed()
	{
		int s, sum = 0;
		for (int i = 1; i < families.size(); i++)
		{
			s = (families[i].score / score) * size;
			sum += s;
			families[i].breed(s);
		}
		// the best familiy gets the extra
		families.front().breed(size - sum);
	}

	void evaluate()
	{
		score = 0;
		for (auto& f : families)
			score += f.evaluate();
		std::sort(families.begin(), families.end(), [](auto const& x, auto const& y){ return x.score > y.score; });
	}
	
	void split()
	{
		Xor_Specie misfits;
		for (auto& f : families)
			for (auto&& c : f.expell_misfits())
				misfits.add(c);
		for (auto& m : misfits)
		{
			for (auto& f : families)
				if (f.is_related_to(m.head()))
				{
					f.add(m);
					goto found_family;
				}
			families.push_back(m);
found_family:;
		}
	}

	void add(Xor_Family const& f)
	{
		families.push_back(f);
	}

	void add(Xor_Creature const& c)
	{
		for (auto& f : families)
			if (f.is_related_to(c))
			{
				f.add(c);
				break;
			}
	}
	
	Xor_Creature& best() 
	{
		return families.front().head();
	}

	void print_results()
	{
		std::cout << "\033[2J\033[H" << std::endl;
		for (auto& f : families)
			f.print_results();
	}

	void print_scores()
	{
		std::cout << "\033[2J\033[H" << std::endl;
		for (auto& f : families)
			f.print_scores();
	}

	Xor_Creature& operator[](int index)
	{
		for (auto& f : families)
			if (index < f.size())
				return f[index];
			else
				index -= f.size();
		return best();
	}

	int _size() const
	{
		int s = 0;
		for (auto& f : families)
			s += f.size();
		return s;
	}
};

SDL_Color black	{ 0,   0,   0, 	 255 };
SDL_Color red	{ 255, 0, 	0, 	 255 };
SDL_Color green	{ 0,   255, 0, 	 255 };
SDL_Color blue	{ 0,   0, 	255, 255 };
SDL_Color white	{ 255, 255,	255, 255 };
SDL_Color yellow{ 255, 255, 0,   255 }; 

int test_mutations()
{
	Xor_Brain brain[3];

	Graphic::Window w{ "mutations test", "font.ttf" };

	// parent 1 controls
	w.add_control(Graphic::Button<[](auto& b){ b.randomly_change_weight(); }, Xor_Brain>	{ "change weight mutation", 	0, 0,   400, 200, red, white, brain[0] });
	w.add_control(Graphic::Button<[](auto& b){ b.randomly_adjust_weight(); }, Xor_Brain>	{ "adjust weight mutation", 	0, 200, 400, 200, red, white, brain[0] });
	w.add_control(Graphic::Button<[](auto& b){ b.randomly_add_connection(); }, Xor_Brain>	{ "add connection mutation", 	0, 400, 400, 200, red, white, brain[0] });
	w.add_control(Graphic::Button<[](auto& b){ b.randomly_split_connection(); }, Xor_Brain>	{ "split connection mutation", 	0, 600, 400, 200, red, white, brain[0] });
	w.add_control(Graphic::Button<[](auto& b){ b.randomly_toggle_connection(); }, Xor_Brain>{ "toggle connection mutation", 0, 800, 400, 200, red, white, brain[0] });
	
	// parent 2 controls
	w.add_control(Graphic::Button<[](auto& b){ b.randomly_change_weight(); }, Xor_Brain>	{ "change weight mutation", 	400, 0,    400, 200, green, white, brain[1] });
	w.add_control(Graphic::Button<[](auto& b){ b.randomly_adjust_weight(); }, Xor_Brain>	{ "adjust weight mutation", 	400, 200,  400, 200, green, white, brain[1] });
	w.add_control(Graphic::Button<[](auto& b){ b.randomly_add_connection(); }, Xor_Brain>	{ "add connection mutation", 	400, 400,  400, 200, green, white, brain[1] });
	w.add_control(Graphic::Button<[](auto& b){ b.randomly_split_connection(); }, Xor_Brain>	{ "split connection mutation", 	400, 600,  400, 200, green, white, brain[1] });
	w.add_control(Graphic::Button<[](auto& b){ b.randomly_toggle_connection(); }, Xor_Brain>{ "toggle connection mutation", 400, 800,  400, 200, green, white, brain[1] });

	// debug controls
	w.add_control(Graphic::Button<[](auto& b){ b.print_genes(); }, Xor_Brain>{ "print genes", 0, 1000, 400, 200, red, white, brain[0] });
	w.add_control(Graphic::Button<[](auto& b){ b.print_nodes(); }, Xor_Brain>{ "print nodes", 0, 1200, 400, 200, red, white, brain[0] });
	
	w.add_control(Graphic::Button<[](auto& b){ b.print_genes(); }, Xor_Brain>{ "print genes", 400, 1000, 400, 200, green, white, brain[1] });
	w.add_control(Graphic::Button<[](auto& b){ b.print_nodes(); }, Xor_Brain>{ "print nodes", 400, 1200, 400, 200, green, white, brain[1] });

	// kid controls
	w.add_control(Graphic::Button<[](auto& p, auto& q, auto& k){ k = p.crossover_with(q); }, Xor_Brain, Xor_Brain, Xor_Brain>{ "crossover", 800, 0, 400, 200, blue, white, brain[0], brain[1], brain[2] });
	
	// brains
	w.add_control(Graphic::Network<[](auto& b) -> auto& { return b; }, Xor_Brain>{ 1200, 0,    1600, 600, red,   yellow, green, red, white, brain[0] });
	w.add_control(Graphic::Network<[](auto& b) -> auto& { return b; }, Xor_Brain>{ 1200, 600,  1600, 600, green, yellow, green, red, white, brain[1] });
	w.add_control(Graphic::Network<[](auto& b) -> auto& { return b; }, Xor_Brain>{ 1200, 1200, 1600, 600, blue,  yellow, green, red, white, brain[2] });

	return w.run();
}

int test_neat()
{
	Xor_Specie xors{ 150 };
	
	Graphic::Window w{ "neat test", "font.ttf" };

	// evolution controls
	w.add_control(Graphic::Button<[](auto& x){ x.evaluate(); }, Xor_Specie>         { "evaluate", 0, 0,   400, 200, red, white, xors });
	w.add_control(Graphic::Button<[](auto& x){ x.split(); x.breed(); }, Xor_Specie>{ "breed", 	 0, 200, 400, 200, red, white, xors });
	w.add_control(Graphic::Button<[](auto& x){ x.mutate(); }, Xor_Specie>          { "mutate",	 0, 400, 400, 200, red, white, xors });
	w.add_control(Graphic::Button<[](auto& x)
	{ 
		x.evaluate();
		for (int i = 0; i < 100; i++)
		{
			x.split();
			x.breed();
			x.mutate(); 
			x.evaluate();
		}
	}, Xor_Specie>          { "evolve for 100 generations",	 0, 600, 400, 200, red, white, xors });
	
	// debug controls
	w.add_control(Graphic::Button<[](auto& x){ x.print_results(); }, Xor_Specie>{ "print results", 0, 800, 400, 200, red, white, xors });
	w.add_control(Graphic::Button<[](auto& x){ x.print_scores(); }, Xor_Specie>{ "print scores", 0, 1000, 400, 200, red, white, xors });

	// brains
	w.add_control(Graphic::Network<[](auto& x) -> auto& { return x[0].brain; }, Xor_Specie>  { 400, 0,   1600, 600, green, yellow, green, red, white, xors });
	w.add_control(Graphic::Network<[](auto& x) -> auto& { return x[1].brain; }, Xor_Specie>{ 400, 600, 1600, 600, blue,  yellow, green, red, white, xors });

	return w.run();
}

int main()
{
//	return test_mutations();
	return test_neat();
}	
