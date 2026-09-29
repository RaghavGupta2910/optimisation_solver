# ============================================================
# Demo applications
# ============================================================

# Common demo utilities
add_library(demo_common
    common/demo_utils.cpp
    common/data_utils.cpp
    common/reporting.cpp
    common/qp3_loader.cpp
)

target_link_libraries(demo_common
    PUBLIC
        mps_reader
        solver_orchestrator
        pdlp_adapter
        qp_engine
)

target_include_directories(demo_common
    PUBLIC
        ${CMAKE_SOURCE_DIR}
)

# ------------------------------------------------------------
# MILP - Supply Chain
# ------------------------------------------------------------
add_executable(supply_chain_milp
    milp/supply_chain/main.cpp
)

target_link_libraries(supply_chain_milp
    PRIVATE
        demo_common
        milp_engine
)

# ------------------------------------------------------------
# PDLP - Transportation
# ------------------------------------------------------------
add_executable(transportation_pdlp
    pdlp/transportation/main.cpp
)

target_link_libraries(transportation_pdlp
    PRIVATE
        demo_common
        pdlp_engine
)

# ------------------------------------------------------------
# QP - Portfolio
# ------------------------------------------------------------
add_executable(portfolio_qp
    qp/portfolio/main.cpp
)

target_link_libraries(portfolio_qp
    PRIVATE
        demo_common
        qp_engine
)

# ------------------------------------------------------------
# SuperADMM - Portfolio
# ------------------------------------------------------------
add_executable(portfolio_superadmm
    qp_superadmm/portfolio/main.cpp
    qp_superadmm/portfolio/convergence.cpp
)

target_link_libraries(portfolio_superadmm
    PRIVATE
        demo_common
        qp_engine
)

# ------------------------------------------------------------
# MIQP - Unit Commitment
# ------------------------------------------------------------
add_executable(unit_commitment_miqp
    miqp/unit_commitment/main.cpp
)

target_link_libraries(unit_commitment_miqp
    PRIVATE
        demo_common
        milp_engine
        qp_engine
)

# ------------------------------------------------------------
# NLP - Reactor Optimization
# ------------------------------------------------------------
add_executable(reactor_optimization_nlp
    nlp/reactor_optimization/main.cpp
)

target_link_libraries(reactor_optimization_nlp
    PRIVATE
        demo_common
)

# ------------------------------------------------------------
# Organize demo targets
# ------------------------------------------------------------
set_target_properties(
    supply_chain_milp
    transportation_pdlp
    portfolio_qp
    portfolio_superadmm
    unit_commitment_miqp
    reactor_optimization_nlp
    PROPERTIES
        FOLDER "Demos"
)

