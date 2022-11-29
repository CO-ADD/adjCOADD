const resizer = document.querySelector(".resizer");
const sidebar = document.querySelector(".sidebar");

resizer.addEventListener("mousedown", (event) => {
    console.log("enve")
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

resizer.addEventListener("dblclick", () => {
    console.log("double clicked")
    sidebar.classList.toggle("not-visible");
    resizer.classList.toggle("toLeft")
})