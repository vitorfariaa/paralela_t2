#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <cstdint>
#include <chrono>
#include <cstring>

#ifdef USE_MPI
#include <mpi.h>
#endif

#define standard_input  std::cin
#define standard_output std::cout

using namespace std;

using Boolean = bool ;
using Size    = std::size_t ;
using String  = std::string ;

using InStream  = std::istream ;
using OutStream = std::ostream ;

template <typename T, typename U>
using Pair = std::pair <T, U> ;

template <typename T, typename C = std::less <T>>
using Set = std::set <T> ;

template <typename T>
using SizeType = typename T :: size_type ;

template <typename C> inline auto
size (const C& x) -> SizeType <C>
{
    return x.size () ;
}

template <typename C> inline auto
at_least_two_elements_in (const C& c) -> Boolean 
{ 
    return size (c) > SizeType <C> (1); 
}

template <typename T> inline auto
first_element (const Set <T>& x) -> T 
{ 
    return *(x.begin ()); 
}

template <typename T> inline auto
second_element (const Set <T>& x) -> T 
{ 
    return *(std::next (x.begin ()));
 }

template <typename T> inline auto
remove (Set <T>& x, const T& e) -> Set <T>& 
{ 
    x.erase (e);
    return x; 
}

template <typename T> inline auto
push (Set <T>& x, const T& e) -> Set <T>& 
{ 
    x.insert (e); 
    return x; 
}

template <typename C> inline auto
empty (const C& x) -> Boolean 
{ 
    return size (x) == SizeType <C> (0); 
}

static double tempo_paralelo = 0.0; // acumula o tempo da parte paralelizavel

// parte sequencial 
// definicoes do enunciado: prefixo, sufixo, overlap, etc

inline auto 
Boolean is_prefix (const String& a, const String& b)
{
    if (size (a) > size (b)) return false;
    if (! (std::mismatch(a.begin(), a.end(), b.begin()).first == a.end())) return false;
    return true;
}

inline auto suffix_from_position (const String& x, SizeType <String> i) -> String 
{ 
    return x.substr (i); 
}

inline auto remove_prefix (const String& x, SizeType <String> n) -> String
{
    if (size (x) > n) return suffix_from_position (x, n);
    return x;
}

auto all_suffixes (const String& x) -> Set <String>
{
    Set <String> ss;
    SizeType <String> n = size (x);
    while (--n) { ss.insert (x.substr (n)); }
    return ss;
}

auto commom_suffix_and_prefix (const String& a, const String& b) -> String
{
    if (empty (a) || empty (b)) return "";
    String x = "";
    for (const String& s : all_suffixes (a)) {
        if (is_prefix (s, b) && size (s) > size (x)) x = s;
    }
    return x;
}

inline auto overlap_value (const String& s, const String& t) -> SizeType <String>
{
    return size (commom_suffix_and_prefix (s, t));
}

auto overlap (const String& s, const String& t) -> String
{
    String c = commom_suffix_and_prefix (s, t);
    return s + remove_prefix (t, size (c));
}

inline auto pop_two_elements_and_push_overlap (Set <String>& ss, const Pair <String, String>& p) -> Set <String>&
{
    ss = remove (ss, p.first);
    ss = remove (ss, p.second);
    ss = push   (ss, overlap (p.first, p.second));
    return ss;
}

auto all_distinct_pairs (const Set <String>& ss) -> Set <Pair <String, String>>
{
    Set <Pair <String, String>> x;
    for (const String& s1 : ss) {
        for (const String& s2 : ss) {
            if (s1 != s2) x.insert (std::make_pair (s1, s2));
        }
    }
    return x;
}

auto highest_overlap_value (const Set <Pair <String, String>>& sp) -> Pair <String, String>
{
    Pair <String, String> x = first_element (sp);
    for (const Pair <String, String>& p : sp) {
        if (overlap_value (p.first, p.second) > overlap_value (x.first, x.second)) {
            x = p;
        }
    }
    return x;
}

// novas funcoes -> codigo paralelizado
// ja que o set deixa os pares ordenados por default, funcao de comparador serve pra desempatar lexicograficamente
static inline auto comparador (const Pair <String,String>& a, const Pair <String,String>& b) -> Boolean
{
    return (a.first < b.first) || (a.first == b.first && a.second < b.second);
}

// pequena melhoria: remove cadeias que sao substrings de outras
static auto remove_redundant_substrings (const Set <String>& ss) -> Set <String>
{
    vector <String> v (ss.begin (), ss.end ());
    Size const n = v.size ();
    vector <Boolean> remove (n, false);

    for (Size i = 0 ; i < n ; ++i) {
        if (remove [i]) continue;
        for (Size j = 0 ; j < n ; ++j) {
            if (i == j || remove [j]) continue;

            if (size (v [i]) >= size (v [j])) {
                if (v [i].find (v [j]) != String :: npos) {
                    remove [j] = true;
                }
            }
        }
    }

    Set <String> out ;
    for (Size i = 0 ; i < n ; ++i) {
        if (!remove [i]) {
            out.insert (v [i]);
        }
    }
    return out;
}

// versao sequencial do resolvedor (usada quando o codigo e compilado sem USE_MPI)
static auto resolvedor_ssp_sequencial (const vector <String>& v) -> Pair <String,String>
{
    auto t0 = chrono :: high_resolution_clock :: now ();

    int const n = (int) v.size ();
    Pair <String,String> melhor = make_pair (v [0], v [1]);
    auto melhor_ov = overlap_value (melhor.first, melhor.second);

    for (int i = 0 ; i < n ; ++i) {
        for (int j = 0 ; j < n ; ++j) {
            if (i == j) continue;

            Pair <String,String> atual = make_pair (v [i], v [j]);
            auto ov = overlap_value (atual.first, atual.second);

            if ( ov > melhor_ov
              || (ov == melhor_ov && comparador (atual, melhor)) )
            {
                melhor    = atual;
                melhor_ov = ov;
            }
        }
    }

    auto t1 = chrono :: high_resolution_clock :: now ();
    tempo_paralelo += chrono :: duration <double> (t1 - t0).count ();
    return melhor;
}

#ifdef USE_MPI

static int mpi_rank = 0 ;
static int mpi_size = 1 ;

// faz o broadcast de um vetor de strings a partir do rank 0
static auto bcast_strings (vector <String>& v) -> void
{
    int root = 0 ;
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);

    int n = (int) v.size ();
    MPI_Bcast (&n, 1, MPI_INT, root, MPI_COMM_WORLD);

    if (mpi_rank != root) {
        v.resize (n);
    }

    if (n == 0) return ;

    vector <int> lens (n);
    if (mpi_rank == root) {
        for (int i = 0 ; i < n ; ++i) {
            lens [i] = (int) v [i].size ();
        }
    }

    MPI_Bcast (lens.data (), n, MPI_INT, root, MPI_COMM_WORLD);

    int total = 0 ;
    for (int i = 0 ; i < n ; ++i) total += lens [i];

    vector <char> buf ((Size) total);
    if (mpi_rank == root) {
        int off = 0 ;
        for (int i = 0 ; i < n ; ++i) {
            std :: memcpy (buf.data () + off, v [i].data (), (Size) lens [i]);
            off += lens [i];
        }
    }

    MPI_Bcast (buf.data (), total, MPI_CHAR, root, MPI_COMM_WORLD);

    if (mpi_rank != root) {
        int off = 0 ;
        for (int i = 0 ; i < n ; ++i) {
            v [i].assign (buf.data () + off, (Size) lens [i]);
            off += lens [i];
        }
    }
}

// mapeia um indice linear k no par (i,j) com i != j
static inline auto linear_to_pair (long long k, int n, int& i, int& j) -> void
{
    i = (int) (k / (n - 1));
    int r = (int) (k % (n - 1));
    j = (r < i) ? r : (r + 1);
}

// resolvedor paralelo usando MPI: divide o espaco de pares (i,j) entre os processos
static auto resolvedor_ssp_mpi (const vector <String>& v) -> Pair <String,String>
{
    int n = (int) v.size ();
    if (n < 2) {
        return make_pair (v [0], v [0]);
    }

    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);

    long long total = 1LL * n * (n - 1);

    long long ini = (total * mpi_rank) / mpi_size;
    long long fim = (total * (mpi_rank + 1)) / mpi_size;

    Pair <String,String> melhor_local ;
    int melhor_local_ov = -1 ;

    auto t0 = chrono :: high_resolution_clock :: now ();

    for (long long k = ini ; k < fim ; ++k) {
        int i, j;
        linear_to_pair (k, n, i, j);

        Pair <String,String> atual = make_pair (v [i], v [j]);
        int ov = (int) overlap_value (atual.first, atual.second);

        if (melhor_local_ov < 0
         || ov > melhor_local_ov
         || (ov == melhor_local_ov && comparador (atual, melhor_local)))
        {
            melhor_local    = atual;
            melhor_local_ov = ov;
        }
    }

    auto t1 = chrono :: high_resolution_clock :: now ();
    double local_time = chrono :: duration <double> (t1 - t0).count ();

    double max_time = 0.0 ;
    MPI_Reduce (&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        tempo_paralelo += max_time;
    }

    struct SimpleBest {
        int ov ;
        int i  ;
        int j  ;
    };

    SimpleBest local_best ;
    local_best.ov = melhor_local_ov;
    local_best.i  = -1;
    local_best.j  = -1;

    if (melhor_local_ov >= 0) {
        // como guardamos apenas indices durante a busca, ja sabemos (i,j)
        // correspondentes ao melhor par local
        int n2 = n;
        for (int i = 0 ; i < n2 ; ++i) {
            for (int j = 0 ; j < n2 ; ++j) {
                if (i == j) continue;
                if (v [i] == melhor_local.first && v [j] == melhor_local.second) {
                    local_best.i = i;
                    local_best.j = j;
                    i = n2; // quebra os dois loops
                    break;
                }
            }
        }
    }

    vector <SimpleBest> all;
    void* recvbuf = nullptr;
    if (mpi_rank == 0) {
        all.resize (mpi_size);
        recvbuf = static_cast <void*> (all.data ());
    }

    MPI_Gather (&local_best, sizeof (SimpleBest), MPI_BYTE,
                recvbuf,    sizeof (SimpleBest), MPI_BYTE,
                0, MPI_COMM_WORLD);

    int best_i = 0, best_j = 1;
    int best_ov = -1;

    if (mpi_rank == 0) {
        for (int r = 0 ; r < mpi_size ; ++r) {
            if (all [r].ov < 0) continue;

            Pair <String,String> atual = make_pair (v [all [r].i], v [all [r].j]);
            int ov = all [r].ov;

            if ( best_ov < 0
              || ov > best_ov
              || (ov == best_ov && comparador (atual, make_pair (v [best_i], v [best_j]))) )
            {
                best_ov = ov;
                best_i  = all [r].i;
                best_j  = all [r].j;
            }
        }
    }

    MPI_Bcast (&best_i, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&best_j, 1, MPI_INT, 0, MPI_COMM_WORLD);

    return make_pair (v [best_i], v [best_j]);
}

#endif // USE_MPI

// funcao de alto nivel usada pelo kernel guloso
static auto resolvedor_ssp (const vector <String>& v) -> Pair <String,String>
{
#ifdef USE_MPI
    return resolvedor_ssp_mpi (v);
#else
    return resolvedor_ssp_sequencial (v);
#endif
}

// faz o equivalente do "return highest_overlap_value (all_distinct_pairs (ss))";
auto pair_of_strings_with_highest_overlap_value (const Set <String>& ss) -> Pair <String, String>
{
    vector <String> v (ss.begin (), ss.end ());
    return resolvedor_ssp (v);
}

// parte sequencial -> nao paralelizada
auto shortest_superstring (Set <String> t) -> String
{
    if (empty (t)) return "";

    // guloso com dependencia entre iteracoes -> sequencial
    while (at_least_two_elements_in (t)) {
        t = pop_two_elements_and_push_overlap(t, pair_of_strings_with_highest_overlap_value (t));
    }
    return first_element (t);
}

inline auto write_string_and_break_line (OutStream& out, String s) -> void 
{ 
    out << s << std::endl; 
}

inline auto read_size (InStream& in) -> Size 
{ 
    Size n; 
    in >> n; 
    return n; 
}

inline auto read_string (InStream& in) -> String { 
    String s; 
    in >> s; 
    return s; 
}

auto read_strings_from_standard_input () -> Set <String>
{
    Set <String> ss;
    SizeType <String> n = read_size (standard_input);
    while (n--) ss.insert (read_string (standard_input));
    return ss;
}

inline auto write_string_to_standard_ouput (const String& s) -> void 
{ 
    write_string_and_break_line (standard_output, s); 
}

auto main (int argc, char const* argv[]) -> int
{
#ifdef USE_MPI
    MPI_Init (&argc, (char***) &argv);

    int rank = 0 ;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    // leitura apenas no rank 0
    vector <String> inicial ;
    if (rank == 0) {
        Set <String> ss_in = read_strings_from_standard_input ();
        ss_in = remove_redundant_substrings (ss_in);
        inicial.assign (ss_in.begin (), ss_in.end ());
    }

    // distribui as strings para todos os processos
    bcast_strings (inicial);

    Set <String> ss ;
    for (const auto& s : inicial) {
        ss.insert (s);
    }

    auto start = chrono :: high_resolution_clock :: now ();
    String resultado = shortest_superstring (ss);
    auto end   = chrono :: high_resolution_clock :: now ();

    double total = chrono :: duration <double> (end - start).count ();

    if (rank == 0) {
        write_string_to_standard_ouput (resultado);
        cout << ((total - tempo_paralelo) / total) * 100 << "%\n";
    }

    MPI_Finalize ();
    return 0;
#else
    auto start = chrono :: high_resolution_clock :: now ();
    Set <String> ss = read_strings_from_standard_input ();
    ss = remove_redundant_substrings (ss);
    write_string_to_standard_ouput (shortest_superstring (ss));
    auto end = chrono :: high_resolution_clock :: now ();

    double total = chrono :: duration <double> (end - start).count ();
    cout << ((total - tempo_paralelo) / total) * 100 << "%\n";
    return 0;
#endif
}
