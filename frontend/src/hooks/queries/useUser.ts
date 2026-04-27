/**
 * TanStack Query hooks for the current User profile/settings/tokens.
 */
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'

import { queryKeys } from '@/lib/queryClient'
import userService from '@/services/user.service'

export function useCurrentUserProfile() {
  return useQuery({
    queryKey: queryKeys.auth.me,
    queryFn: () => userService.getProfile(),
    staleTime: 60_000,
  })
}

export function useUserStats() {
  return useQuery({
    queryKey: ['user', 'stats'],
    queryFn: () => userService.getUserStats(),
    staleTime: 60_000,
  })
}

export function useUserActivity(limit = 10) {
  return useQuery({
    queryKey: ['user', 'activity', limit],
    queryFn: () => userService.getUserActivity(limit),
    staleTime: 30_000,
  })
}

export function useUserSettings() {
  return useQuery({
    queryKey: ['user', 'settings'],
    queryFn: () => userService.getSettings(),
    staleTime: 5 * 60_000,
  })
}

export function useUserNotificationSettings() {
  return useQuery({
    queryKey: ['user', 'notification-settings'],
    queryFn: () => userService.getNotificationSettings(),
    staleTime: 5 * 60_000,
  })
}

export function useApiTokens() {
  return useQuery({
    queryKey: ['user', 'api-tokens'],
    queryFn: () => userService.getApiTokens(),
  })
}

// ----- mutations -----

export function useUpdateProfile() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (data: any) => userService.updateProfile(data),
    onSuccess: (data) => {
      queryClient.setQueryData(queryKeys.auth.me, data)
    },
  })
}

export function useChangePassword() {
  return useMutation({
    mutationFn: (data: any) => userService.changePassword(data),
  })
}

export function useUploadAvatar() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (file: File) => userService.uploadAvatar(file),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.auth.me })
    },
  })
}

export function useUpdateUserSettings() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (data: any) => userService.updateSettings(data),
    onSuccess: (data) => {
      queryClient.setQueryData(['user', 'settings'], data)
    },
  })
}

export function useUpdateUserNotificationSettings() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (data: any) => userService.updateNotificationSettings(data),
    onSuccess: (data) => {
      queryClient.setQueryData(['user', 'notification-settings'], data)
    },
  })
}

export function useCreateApiToken() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (data: { name: string; description?: string }) => userService.createApiToken(data),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['user', 'api-tokens'] })
    },
  })
}

export function useDeleteApiToken() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (tokenId: string) => userService.deleteApiToken(tokenId),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['user', 'api-tokens'] })
    },
  })
}
