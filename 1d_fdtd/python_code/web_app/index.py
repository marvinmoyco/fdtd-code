import streamlit as st


pages = {
    "Index": [
        st.Page("home.py", title="Home"),
        st.Page("add_sim.py", title="Add new simulation"),
    ],

}

pg = st.navigation(pages)
pg.run()
