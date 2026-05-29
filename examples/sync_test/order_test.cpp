// =====================================================================
//  order_test.cpp  --  direct proof that the kernel preserves source
//  ordering within update(), as required by Section 4.2.2.
//
//  The fundamental claim of 4.2.2 is: the modeler decides whether
//  discrete operations execute before or after continuous derivatives
//  by where they place State::sample() blocks in the source.  The
//  kernel must respect that ordering exactly -- it must NOT reorder
//  the code, batch sample() blocks separately, or otherwise interfere.
//
//  This test instruments update() with a log buffer and confirms that
//  the line ordering observed at runtime matches the source ordering.
// =====================================================================

#include "../../osk/osk.h"
#include <cstdio>
#include <string>
#include <vector>

using namespace osk;

class OrderModel : public Block {
public:
    double x, xd;
    std::vector<std::string> log;

    OrderModel() {
        x = 0.0; xd = 0.0;
        addIntegrator(x, xd);
    }
    void update() override {
        // Only record on the boundary where sample() can fire.
        if (State::stepstart) {
            log.push_back("A: continuous before-sample");
        }
        if (State::sample(0.1)) {
            log.push_back("B: sample block 1");
        }
        if (State::stepstart) {
            log.push_back("C: continuous mid");
        }
        if (State::sample(0.5)) {
            log.push_back("D: sample block 2");
        }
        if (State::stepstart) {
            log.push_back("E: continuous after-sample");
        }
        xd = 0.0;
    }
};

int main() {
    OrderModel* m = new OrderModel();
    std::vector<Block*> stage = { m };
    std::vector< std::vector<Block*> > stages = { stage };
    double dts[] = { 0.01 };
    Sim sim(dts, 0.5, stages);
    sim.run();

    // Expected ordering at every boundary:
    //   t = 0.0   : both samples fire    -> A, B, C, D, E
    //   t = 0.1   : only 0.1-Hz fires    -> A, B, C, E
    //   t = 0.2   : only 0.1-Hz fires    -> A, B, C, E
    //   t = 0.3   : only 0.1-Hz fires    -> A, B, C, E
    //   t = 0.4   : only 0.1-Hz fires    -> A, B, C, E
    //   t = 0.5   : both fire            -> A, B, C, D, E  (also twice
    //                                       because of ticklast duplicate)
    //
    // What MUST be true regardless of the precise stepping is:
    //   - within any contiguous run of log entries between two A's,
    //     the entries appear in source order: A < B? < C < D? < E
    //   - if a sample line appears, it's where the source puts it
    //     (B between A and C; D between C and E).

    int errors = 0;
    int boundaries = 0;
    std::size_t i = 0;
    while (i < m->log.size()) {
        if (m->log[i] != "A: continuous before-sample") {
            std::printf("FAIL: boundary should start with A, got '%s'\n",
                        m->log[i].c_str());
            errors++;
            break;
        }
        // Walk this boundary's entries up to the next A (or end).
        std::vector<std::string> got;
        got.push_back(m->log[i++]);
        while (i < m->log.size() && m->log[i] != "A: continuous before-sample") {
            got.push_back(m->log[i++]);
        }
        // Expected slot order: A, [B?], C, [D?], E.
        std::vector<std::string> expected_skeleton = {
            "A: continuous before-sample",
            "B: sample block 1",
            "C: continuous mid",
            "D: sample block 2",
            "E: continuous after-sample"
        };
        // Walk `got` against the skeleton: B and D are optional, the
        // others must appear and must appear in this order.
        std::size_t s = 0;
        for (const auto& line : got) {
            // advance skeleton until we find a match
            bool matched = false;
            while (s < expected_skeleton.size()) {
                if (line == expected_skeleton[s]) {
                    matched = true;
                    s++;
                    break;
                }
                // skip an optional slot only if it's B or D
                if (expected_skeleton[s] == "B: sample block 1" ||
                    expected_skeleton[s] == "D: sample block 2") {
                    s++;
                } else {
                    break;
                }
            }
            if (!matched) {
                std::printf("FAIL: out-of-order entry '%s'\n", line.c_str());
                errors++;
                break;
            }
        }
        boundaries++;
    }

    std::printf("Inspected %d boundaries, %d errors.\n", boundaries, errors);
    std::printf("Total log entries: %zu\n", m->log.size());
    std::printf("\nFirst 12 entries (t=0 + t=0.1 boundaries):\n");
    for (std::size_t k = 0; k < m->log.size() && k < 12; ++k) {
        std::printf("  [%2zu] %s\n", k, m->log[k].c_str());
    }
    delete m;
    return errors == 0 ? 0 : 1;
}
