$(document).ready(function () {
    $(document).on("click", ".edit_detailtable", function () {
        console.log("loading editableTable...")
        var selectbox = $(".select")
        var detail = $(".detail")
        selectbox.addClass("visible")
        $('.updateformshow').addClass("visible")
        detail.addClass("not-visible")
        $('.editablechoices').removeClass("editablechoices")
        $('.editable').removeClass("editable")

    });

    // $(document).keypress(
    //     function (event) {
    //         if (event.which == '13') {
    //             event.preventDefault();
    //         }
    //     });
})

