const hidefilter = document.getElementById("hidefilter")


hidefilter.addEventListener("click", () => {
    hidefilter.classList.toggle("text-danger")
    sidebar.classList.toggle("not-visible");
    resizer.classList.toggle("toLeft")
})