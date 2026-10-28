from rest_framework.permissions import BasePermission


class IsOwner(BasePermission):
    def has_object_permission(self, request, view, obj):
        return obj.user == request.user


class IsAppUser(BasePermission):
    """
    Allow access only to authenticated users with application access
    (ApplicationUser.is_appuser is True).
    """
    def has_permission(self, request, view):
        user = request.user
        return bool(
            user
            and user.is_authenticated
            and getattr(user, 'is_appuser', False)
        )

