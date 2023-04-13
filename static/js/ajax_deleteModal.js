function ajaxDelete(deleteUrl, csrfToken, successCallback, errorCallback) {
    $.ajax({
        url: deleteUrl,
        type: 'POST',
        data: {
            csrfmiddlewaretoken: csrfToken
        },
        success: function (data) {
            if (typeof successCallback === 'function') {
                successCallback(data);
            }
        },
        error: function (xhr, errmsg, err) {
            if (typeof errorCallback === 'function') {
                errorCallback(xhr, errmsg, err);
            } else {
                console.error(xhr.status + ": " + xhr.responseText);
            }
        }
    });
}