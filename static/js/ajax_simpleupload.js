function ajaxSimpleUpload(inputSelector, inputfield, url, fileNameInputSelector, fileTypeInputSelector) {

    const csrftoken = getCookie("csrftoken");
    $(inputSelector).on('change', function () {
        console.log(url)
        var fileInput = $(inputSelector)[0];
        var file = fileInput.files[0];
        var formData = new FormData();
        formData.append(inputfield, file);
        $.ajax({
            url: url,
            type: 'POST',
            data: formData,
            processData: false,
            contentType: false,
            headers: { "X-CSRFToken": csrftoken},
            success: function (data) {
                console.log("success")
                $(fileNameInputSelector).val(data.name);
                console.log(data.name)
                console.log( $(fileNameInputSelector))
                $(fileTypeInputSelector).val(data.type);
            },
            error: function (xhr, errmsg, err) {
                if (typeof errorCallback === 'function') {
                    errorCallback(xhr, errmsg, err);
                } else {
                    console.error(xhr.status + ": " + xhr.responseText);
                }
            }
        });
    });
}