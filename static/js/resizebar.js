const resizer = document.querySelector(".resizer");
const sidebar = document.querySelector(".sidebar");
const hidefilter = document.getElementById("hidefilter")


resizer.addEventListener("mousedown", (event) => {
    console.log(hidefilter)
    document.addEventListener("mousemove", resize, false);
    document.addEventListener(
        "mouseup",
        () => {
            document.removeEventListener("mousemove", resize, false);
        },
        false
    );
});

function resize(e) {
    var sdsize = Math.min(500, e.x);
    var size = `${sdsize}px`;
    sidebar.style.flexBasis = size;
}

hidefilter.addEventListener("click", () => {
    hidefilter.classList.toggle("text-danger")
    sidebar.classList.toggle("not-visible");
    resizer.classList.toggle("toLeft")
})